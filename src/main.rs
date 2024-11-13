use rust_htslib::bam::{self, Read, Record};
use std::fs::File;
use std::io::{self, BufRead};
use intervaltree::IntervalTree;
use std::collections::HashMap;
use rayon::prelude::*;
use crossbeam::channel::{bounded, Receiver, Sender};
use clap::{Arg, Command};
use std::error::Error;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Arc;
use std::cmp::Ordering as CmpOrdering;

#[derive(Debug)]
enum ProcessError {
    IoError(io::Error),
    BamError(rust_htslib::errors::Error),
    SendError(String),
}

impl From<io::Error> for ProcessError {
    fn from(err: io::Error) -> ProcessError {
        ProcessError::IoError(err)
    }
}

impl From<rust_htslib::errors::Error> for ProcessError {
    fn from(err: rust_htslib::errors::Error) -> ProcessError {
        ProcessError::BamError(err)
    }
}

impl std::error::Error for ProcessError {}

impl std::fmt::Display for ProcessError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ProcessError::IoError(e) => write!(f, "IO error: {}", e),
            ProcessError::BamError(e) => write!(f, "BAM error: {}", e),
            ProcessError::SendError(e) => write!(f, "Send error: {}", e),
        }
    }
}

#[derive(Clone)]
struct OrderedRecord {
    index: usize,
    record: Record,
}

fn read_bed_file(bed_file: &str) -> Result<HashMap<String, IntervalTree<u64, ()>>, Box<dyn Error>> {
    let mut interval_data: HashMap<String, Vec<(std::ops::Range<u64>, ())>> = HashMap::new();
    let file = File::open(bed_file)?;
    
    // 首先收集所有区间
    for line in io::BufReader::new(file).lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 3 {
            let chrom = fields[0].to_string();
            let start: u64 = fields[1].parse()?;
            let end: u64 = fields[2].parse()?;
            interval_data.entry(chrom)
                .or_default()
                .push((start..end, ()));
        }
    }
    
    // 然后为每个染色体创建 IntervalTree
    let intervals_map = interval_data.into_iter()
        .map(|(chrom, intervals)| {
            (chrom, IntervalTree::from_iter(intervals))
        })
        .collect();
    
    Ok(intervals_map)
}

fn process_records(
    intervals_map: &HashMap<String, IntervalTree<u64, ()>>,
    tid_to_name: HashMap<i32, String>,
    receiver: Receiver<OrderedRecord>,
    sender: Sender<OrderedRecord>,
    processed_count: Arc<AtomicUsize>,
) {
    receiver.into_iter().par_bridge().for_each_with(sender, |s, mut ordered_record| {
        let tid = ordered_record.record.tid();
        if tid >= 0 {
            if let Some(tid_name) = tid_to_name.get(&tid) {
                if let Some(tree) = intervals_map.get(tid_name) {
                    let record_start = ordered_record.record.pos() as u64;
                    let record_end = ordered_record.record.cigar().end_pos() as u64;
                    if tree.query(record_start..record_end).next().is_some() {
                        if ordered_record.record.mapq() < 30 {
                            ordered_record.record.set_mapq(60);
                        }
                    }
                }
            }
        }
        if let Err(e) = s.send(ordered_record) {
            eprintln!("Error sending record: {}", e);
        }
        processed_count.fetch_add(1, Ordering::Relaxed);
    });
}

fn update_bam_mapq(bam_file: &str, bed_file: &str, output_bam_file: &str) -> Result<(), Box<dyn Error>> {
    println!("Reading BED file...");
    let intervals_map = read_bed_file(bed_file)?;
    println!("Found {} chromosomes in BED file", intervals_map.len());

    println!("Opening BAM file...");
    let mut bam = bam::Reader::from_path(bam_file)?;
    let header = bam::Header::from_template(bam.header());
    let bam_header = bam.header();
    
    println!("Creating chromosome name mapping...");
    let mut tid_to_name = HashMap::new();
    for tid in 0..bam_header.target_count() {
        if let Ok(name) = std::str::from_utf8(bam_header.tid2name(tid)) {
            tid_to_name.insert(tid as i32, name.to_string());
        }
    }
    println!("Found {} chromosomes in BAM file", tid_to_name.len());
    
    println!("Creating output BAM file...");
    let mut output_bam = bam::Writer::from_path(output_bam_file, &header, bam::Format::Bam)?;
    output_bam.set_threads(4)?;

    let (record_sender, record_receiver) = bounded(10000);
    let (result_sender, result_receiver) = bounded(10000);
    
    let processed_count = Arc::new(AtomicUsize::new(0));
    let processed_count_clone = processed_count.clone();

    println!("Starting processing...");
    crossbeam::scope(|scope| -> Result<(), Box<dyn Error>> {
        // 处理记录的线程
        scope.spawn(|_| {
            process_records(
                &intervals_map,
                tid_to_name,
                record_receiver,
                result_sender,
                processed_count_clone,
            );
        });

        // 读取记录的线程
        let read_handle = scope.spawn(|_| -> Result<usize, ProcessError> {
            let mut count = 0;
            let mut records = Vec::with_capacity(1000);
            
            for result in bam.records() {
                match result {
                    Ok(record) => {
                        records.push(OrderedRecord {
                            index: count,
                            record,
                        });
                        count += 1;
                        
                        if records.len() >= 1000 {
                            for record in records.drain(..) {
                                if let Err(e) = record_sender.send(record) {
                                    return Err(ProcessError::SendError(e.to_string()));
                                }
                            }
                            if count % 1_000_000 == 0 {
                                println!("Read {} records", count);
                            }
                        }
                    }
                    Err(e) => return Err(ProcessError::BamError(e)),
                }
            }
            
            // 发送剩余的记录
            for record in records.drain(..) {
                if let Err(e) = record_sender.send(record) {
                    return Err(ProcessError::SendError(e.to_string()));
                }
            }
            
            drop(record_sender);
            Ok(count)
        });

        // 写入记录的线程 - 使用排序的批量写入
        let mut write_count = 0;
        let mut next_index = 0;
        let mut pending_records: Vec<OrderedRecord> = Vec::new();

        for ordered_record in result_receiver {
            if ordered_record.index == next_index {
                // 如果是期望的下一条记录，直接写入
                output_bam.write(&ordered_record.record)?;
                write_count += 1;
                next_index += 1;

                // 检查是否有待处理的记录可以写入
                pending_records.sort_by_key(|r: &OrderedRecord| r.index);
                while let Some(pos) = pending_records.iter().position(|r| r.index == next_index) {
                    let record = pending_records.remove(pos);
                    output_bam.write(&record.record)?;
                    write_count += 1;
                    next_index += 1;
                }
            } else {
                // 如果不是期望的下一条记录，加入待处理队列
                pending_records.push(ordered_record);
            }

            if write_count % 1_000_000 == 0 {
                println!("Wrote {} records", write_count);
            }
        }

        // 处理剩余的记录
        pending_records.sort_by_key(|r| r.index);
        for record in pending_records {
            output_bam.write(&record.record)?;
            write_count += 1;
        }

        if let Ok(total_read) = read_handle.join() {
            match total_read {
                Ok(count) => {
                    println!("Total records read: {}", count);
                    println!("Total records processed: {}", processed_count.load(Ordering::Relaxed));
                    println!("Total records written: {}", write_count);
                }
                Err(e) => eprintln!("Error in read thread: {}", e),
            }
        }

        Ok(())
    }).unwrap()?;

    println!("Processing complete!");
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let matches = Command::new("BAM MAPQ Updater")
        .version("1.0")
        .author("Your Name <your.email@example.com>")
        .about("Updates MAPQ values in a BAM file based on regions defined in a BED file")
        .arg(Arg::new("bam")
            .short('b')
            .long("bam")
            .value_name("BAM_FILE")
            .help("Input BAM file")
            .required(true))
        .arg(Arg::new("bed")
            .short('d')
            .long("bed")
            .value_name("BED_FILE")
            .help("Input BED file")
            .required(true))
        .arg(Arg::new("output")
            .short('o')
            .long("output")
            .value_name("OUTPUT_BAM_FILE")
            .help("Output BAM file")
            .required(true))
        .get_matches();

    let bam_file = matches.get_one::<String>("bam").unwrap();
    let bed_file = matches.get_one::<String>("bed").unwrap();
    let output_bam_file = matches.get_one::<String>("output").unwrap();

    update_bam_mapq(bam_file, bed_file, output_bam_file)
}
