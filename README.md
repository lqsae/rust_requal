# BAM MAPQ Updater

一个高性能的工具，用于根据 BED 文件中定义的区域更新 BAM 文件中的 MAPQ 值。

## 安装

克隆仓库
```bash
git clone [repository_url]
```
进入目录
```bash
cd bam-mapq-updater
```

编译
```bash
cargo build --release
```
编译后的程序位于
./target/release/rust_requal

## 使用方法

```bash
./target/release/rust_requal -b input.bam -d target.bed -o output.bam
```

## 参数说明：
- `-b, --bam <BAM_FILE>`: 输入的 BAM 文件
- `-d, --bed <BED_FILE>`: 包含目标区域的 BED 文件
- `-o, --output <OUTPUT_BAM_FILE>`: 输出的 BAM 文件

## 工作原理

1. 读取 BED 文件中定义的基因组区域
2. 使用区间树数据结构高效存储和查询这些区域
3. 并行处理 BAM 文件中的记录
4. 对于落在 BED 文件定义区域内且 MAPQ < 30 的记录，将其 MAPQ 值更新为 60
5. 使用多线程压缩写入输出文件

## 性能优化

- 使用 rayon 进行并行处理
- 批量读取和写入记录（每批 1000 条记录）
- 使用 crossbeam 通道进行高效的线程间通信（缓冲区大小 10000）
- 4 线程 BAM 压缩
- 使用 drain 优化内存使用


## 依赖项

- rust-htslib (0.40): BAM 文件处理
- intervaltree (0.2): 区间查询
- crossbeam (0.8): 线程间通信
- clap (4.3): 命令行参数解析
- rayon (1.7): 并行计算

## 系统要求

- 64 位操作系统
- Rust 1.56.0 或更高版本
- 足够的内存来存储 BED 文件定义的区间
- 推荐多核处理器以获得最佳性能

## 进度报告

程序会在处理过程中报告进度：
- 每读取 100 万条记录报告一次
- 每写入 100 万条记录报告一次
- 完成时报告总读取、处理和写入的记录数

## 错误处理

程序包含完整的错误处理机制：
- 文件读写错误
- BAM 格式错误
- 线程通信错误
- 内存分配错误

## 注意事项

1. 确保有足够的磁盘空间存储输出文件
2. 对于大文件处理，建议使用 SSD 以获得更好的 I/O 性能
3. 程序会自动使用多线程，无需手动配置

## 许可证

MIT License

## 贡献

欢迎提交 Issue 和 Pull Request！

## 更新日志

### v1.0.0
- 初始发布
- 支持基本的 MAPQ 更新功能
- 多线程并行处理
- 进度报告
- 保证输出顺序

