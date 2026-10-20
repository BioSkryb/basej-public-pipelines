// bam-insort: parallel, memory-bounded coordinate sort + BAI index for BAM.
//
// Reads query-grouped BAM from stdin (e.g. the dupblaster output), coordinate
// sorts by (reference_sequence_id, alignment_start) — unmapped reads last — and
// writes a coordinate-sorted BGZF BAM plus its .bai. A Rust replacement for
// `samtools sort | samtools index`.
//
// ALGORITHM
// ---------
// Records are handled as OPAQUE raw BAM record blocks (u32 block_size + bytes):
// we never parse/re-encode them, only peek refID (bytes 0..4) + pos (bytes 4..8)
// to compute the sort key and route them. This makes the byte content of the
// output identical to the input records (safe for every downstream caller).
//
//   * Fast path (fits in --max-mem-gb): buffer all records, one rayon parallel
//     sort, one multithreaded-BGZF write. No temp files.
//
//   * Bucket path (exceeds budget — full-depth WGS): coordinate sort is naturally
//     partitionable by genomic interval, so instead of an external merge we route
//     each record into one of N ordered genomic buckets (by cumulative genome
//     offset), spilling buckets to disk when over budget. At EOF each bucket is
//     sorted + BGZF-compressed INDEPENDENTLY and IN PARALLEL (rayon), then the
//     compressed bucket blobs are byte-CONCATENATED in genomic order (BGZF is
//     concatenable) with a single EOF marker. This removes the single-threaded
//     k-way merge that made the external path a serial bottleneck: the "merge"
//     becomes N parallel sorts + a byte copy.
//
// The BAI is built in a final decompression-bound pass over the written BAM,
// following the canonical noodles pattern (bam/examples/bam_index.rs), so the
// index is compatible with `samtools index` (validated: idxstats + region
// queries identical).
//
// Usage: ... | bam-insort -o out.bam [--threads N] [--max-mem-gb G] [--tmp DIR] [--buckets N]
use std::env;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::num::NonZeroUsize;

use noodles_bam as bam;
use noodles_bam::bai;
use noodles_bgzf as bgzf;
use noodles_csi::binning_index::{index::reference_sequence::bin::Chunk, Indexer};
use noodles_sam::alignment::Record as _;
use noodles_sam::Header;
use rayon::prelude::*;

/// Standard 28-byte BGZF EOF marker (SAM spec). Appended once at the very end of
/// the concatenated output; NOT written between bucket blobs.
const BGZF_EOF: &[u8] = &[
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02,
    0x00, 0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

/// Pack (refID, pos) into an i64 sort key: high 32 bits = refID (unmapped -> i32::MAX
/// so it sorts last), low 32 bits = 0-based pos. Ascending i64 order == (refID, pos).
fn key_of(data: &[u8]) -> i64 {
    let refid = i32::from_le_bytes([data[0], data[1], data[2], data[3]]);
    let pos = i32::from_le_bytes([data[4], data[5], data[6], data[7]]);
    let tid_norm = if refid < 0 { i32::MAX } else { refid };
    ((tid_norm as i64) << 32) | ((pos as u32) as i64)
}

fn main() {
    if let Err(e) = run() {
        eprintln!("bam-insort error: {e}");
        std::process::exit(1);
    }
}

struct Args {
    out_path: String,
    threads: usize,
    max_mem_bytes: u64,
    tmp_dir: String,
    buckets: usize,
}

fn parse_args() -> Args {
    let mut out_path = String::new();
    let mut threads = num_cpus_env();
    let mut max_mem_gb: Option<f64> = None;
    let mut tmp_dir = env::var("TMPDIR").unwrap_or_else(|_| ".".to_string());
    let mut buckets = 512usize;
    let args: Vec<String> = env::args().collect();
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-o" | "--out" => {
                out_path = args[i + 1].clone();
                i += 2;
            }
            "-t" | "--threads" => {
                threads = args[i + 1].parse().unwrap_or(threads);
                i += 2;
            }
            "--max-mem-gb" => {
                max_mem_gb = args[i + 1].parse().ok();
                i += 2;
            }
            "--tmp" => {
                tmp_dir = args[i + 1].clone();
                i += 2;
            }
            "--buckets" => {
                buckets = args[i + 1].parse().unwrap_or(buckets).max(1);
                i += 2;
            }
            "-h" | "--help" => {
                eprintln!(
                    "usage: bam-insort -o <out.bam> [--threads N] [--max-mem-gb G] [--tmp DIR] [--buckets N]\n\
                     reads query-grouped BAM from stdin; writes coordinate-sorted BAM + .bai"
                );
                std::process::exit(0);
            }
            _ => i += 1,
        }
    }
    if out_path.is_empty() {
        eprintln!("usage: bam-insort -o <out.bam> [--threads N] [--max-mem-gb G] [--tmp DIR] [--buckets N]");
        std::process::exit(2);
    }
    let max_mem_bytes = match max_mem_gb {
        Some(g) => (g * 1024.0 * 1024.0 * 1024.0) as u64,
        None => (mem_available_bytes() as f64 * 0.60) as u64,
    };
    Args {
        out_path,
        threads,
        max_mem_bytes,
        tmp_dir,
        buckets,
    }
}

/// One genomic bucket's in-RAM records + optional spill file (raw BAM blocks).
/// The spill file is opened/appended/closed per flush round (see `spill_all`), so
/// we never hold thousands of file descriptors open at once (allows many buckets).
struct Bucket {
    mem: Vec<(i64, Vec<u8>)>,
    spill_path: Option<String>,
}
impl Bucket {
    fn new() -> Self {
        Bucket {
            mem: Vec::new(),
            spill_path: None,
        }
    }
}

/// Genome layout for bucketing: cumulative start offset of each reference and
/// the total genome length, so any (refID,pos) maps to a global genomic offset.
struct GenomeLayout {
    cumulative: Vec<i64>,
    total: i64,
    n: usize, // number of mapped buckets; unmapped bucket = index n
}

impl GenomeLayout {
    fn new(header: &Header, n: usize) -> Self {
        let mut cumulative = Vec::with_capacity(header.reference_sequences().len());
        let mut acc: i64 = 0;
        for (_, rs) in header.reference_sequences() {
            cumulative.push(acc);
            acc += usize::from(rs.length()) as i64;
        }
        GenomeLayout {
            cumulative,
            total: acc.max(1),
            n,
        }
    }

    /// Bucket index for a record's sort key. Mapped -> 0..n-1 by global genomic
    /// offset (monotonic in (refID,pos)); unmapped -> n (written last).
    fn bucket_of(&self, key: i64) -> usize {
        let tid = (key >> 32) as i32;
        if tid == i32::MAX {
            return self.n; // unmapped
        }
        let pos = (key & 0xffff_ffff) as i64;
        let g = self.cumulative[tid as usize] + pos;
        let b = ((g as i128 * self.n as i128) / self.total as i128) as usize;
        b.min(self.n - 1)
    }
}

/// Read one raw BAM alignment record (u32 block_size + block_size bytes) into
/// `buf`. Returns Ok(false) at a clean end of stream.
fn read_raw_record<R: Read>(r: &mut R, buf: &mut Vec<u8>) -> io::Result<bool> {
    let mut len = [0u8; 4];
    match r.read_exact(&mut len) {
        Ok(()) => {}
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(false),
        Err(e) => return Err(e),
    }
    let n = u32::from_le_bytes(len) as usize;
    buf.resize(n, 0);
    r.read_exact(buf)?;
    Ok(true)
}

/// Write one raw record (u32 block_size + bytes) to a BGZF/opaque writer.
fn write_raw_record<W: Write>(w: &mut W, data: &[u8]) -> io::Result<()> {
    w.write_all(&(data.len() as u32).to_le_bytes())?;
    w.write_all(data)
}

fn run() -> io::Result<()> {
    let args = parse_args();
    let workers = NonZeroUsize::new(args.threads.max(1)).unwrap();
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .ok();

    // per-record RAM overhead: (i64 key, Vec<u8> ptr/len/cap) + allocator slack
    let rec_overhead: u64 = 8 + 24 + 16;

    // ---- read header (multithreaded BGZF decode), then read raw records ----
    let stdin = BufReader::new(File::open("/dev/stdin")?);
    let decoder = bgzf::MultithreadedReader::with_worker_count(workers, stdin);
    let mut reader = bam::io::Reader::from(decoder);
    let mut header = reader.read_header()?;
    set_coordinate_sorted(&mut header);
    let layout = GenomeLayout::new(&header, args.buckets);
    let mut raw = reader.into_inner(); // bgzf stream positioned at first record

    const PHASE1_LOG_EVERY: u64 = 50_000_000;
    let t0 = std::time::Instant::now();
    let mut mon = ProcMon::new();
    let mut all: Vec<(i64, Vec<u8>)> = Vec::new(); // fast-path accumulator
    let mut buckets: Option<Vec<Bucket>> = None; // Some(..) once we switch to bucket mode
    let mut buf_bytes: u64 = 0;
    let mut total: u64 = 0;
    let mut spill_rounds: u64 = 0;
    let mut next_log: u64 = PHASE1_LOG_EVERY;
    let mut data = Vec::new();

    while read_raw_record(&mut raw, &mut data)? {
        total += 1;
        let key = key_of(&data);
        let cost = data.len() as u64 + rec_overhead;
        let rec = (key, std::mem::take(&mut data));

        if let Some(bk) = buckets.as_mut() {
            bk[layout.bucket_of(key)].mem.push(rec);
            buf_bytes += cost;
            if buf_bytes >= args.max_mem_bytes {
                spill_all(bk, &args.tmp_dir)?;
                spill_rounds += 1;
                eprintln!(
                    "[bam-insort] phase1 spill #{spill_rounds} @ {}M reads | RSS {:.1} GB | {:.0} cores",
                    total / 1_000_000,
                    rss_gb(),
                    mon.cores()
                );
                buf_bytes = 0;
            }
        } else {
            all.push(rec);
            buf_bytes += cost;
            if buf_bytes >= args.max_mem_bytes {
                // switch to bucket mode: distribute what we have, then spill
                eprintln!(
                    "[bam-insort] budget {:.1} GB reached at {}M reads -> bucket mode ({} buckets) | RSS {:.1} GB",
                    args.max_mem_bytes as f64 / 1e9,
                    total / 1_000_000,
                    args.buckets,
                    rss_gb()
                );
                let mut bk: Vec<Bucket> = (0..=args.buckets).map(|_| Bucket::new()).collect();
                for r in all.drain(..) {
                    bk[layout.bucket_of(r.0)].mem.push(r);
                }
                all = Vec::new(); // free the fast-path accumulator's capacity
                spill_all(&mut bk, &args.tmp_dir)?;
                spill_rounds += 1;
                buckets = Some(bk);
                buf_bytes = 0;
            }
        }

        if total >= next_log {
            let secs = t0.elapsed().as_secs_f64();
            eprintln!(
                "[bam-insort] phase1 read {}M reads in {:.0}s ({:.1}M/s) | RSS {:.1} GB | {:.0} cores | {} | spills {}",
                total / 1_000_000,
                secs,
                total as f64 / 1e6 / secs,
                rss_gb(),
                mon.cores(),
                if buckets.is_some() { "bucket" } else { "in-RAM" },
                spill_rounds
            );
            next_log += PHASE1_LOG_EVERY;
        }
    }
    eprintln!(
        "[bam-insort] phase1 DONE: {total} reads in {:.0}s | RSS {:.1} GB | {} | spill rounds {}",
        t0.elapsed().as_secs_f64(),
        rss_gb(),
        if buckets.is_some() { "bucket mode" } else { "in-RAM (fast path)" },
        spill_rounds
    );

    match buckets {
        None => write_fast_path(&args.out_path, &header, &mut all, workers)?,
        Some(mut bk) => {
            // flush any in-RAM tail so every bucket is fully on disk
            spill_all(&mut bk, &args.tmp_dir)?;
            write_bucket_path(&args.out_path, &header, bk, &args.tmp_dir, workers)?;
        }
    }

    // ---- build BAI (decompression-bound final pass; proven noodles path) ----
    let t3 = std::time::Instant::now();
    write_bai(&args.out_path, workers)?;
    eprintln!(
        "[bam-insort] phase3 indexed {}.bai in {:.1}s",
        args.out_path,
        t3.elapsed().as_secs_f64()
    );
    eprintln!(
        "[bam-insort] ALL DONE: {} reads in {:.0}s total | peak RSS {:.1} GB",
        total,
        t0.elapsed().as_secs_f64(),
        peak_rss_gb()
    );
    Ok(())
}

/// Flush every bucket's in-RAM records to its spill file, freeing RAM. The spill
/// file is opened in append mode and closed within this call, so at most one
/// spill file descriptor is open at a time regardless of the bucket count.
fn spill_all(buckets: &mut [Bucket], tmp_dir: &str) -> io::Result<()> {
    use std::fs::OpenOptions;
    // Spill buckets in PARALLEL: zstd-compressing ~29 GB on a single thread stalled
    // phase-1 intake (which, in the fused pipe, backpressures the aligner). Buckets are
    // independent, so we compress+write them across all cores using the CPU that's
    // otherwise idle while bam-insort is just buffering the aligner's stream. Each rayon
    // task holds at most one spill file open, so FDs stay bounded.
    buckets
        .par_iter_mut()
        .enumerate()
        .try_for_each(|(i, b)| -> io::Result<()> {
            if b.mem.is_empty() {
                return Ok(());
            }
            if b.spill_path.is_none() {
                b.spill_path = Some(format!("{}/bam-insort.bucket.{}.zst", tmp_dir, i));
            }
            let path = b.spill_path.as_ref().unwrap();
            // Append one zstd frame per round (zstd decodes concatenated frames on read).
            // Level 1 shrinks spilled bytes ~5-6x -> far less phase-1 write + phase-2 read I/O.
            let f = OpenOptions::new().create(true).append(true).open(path)?;
            let mut enc = zstd::stream::write::Encoder::new(BufWriter::new(f), 1)?;
            for (_, d) in b.mem.drain(..) {
                write_raw_record(&mut enc, &d)?;
            }
            let mut bw = enc.finish()?; // finalize the zstd frame
            bw.flush()?;
            Ok(())
        })
}

/// Fast path: all records fit in RAM. One parallel sort + one multithreaded
/// BGZF write (records copied as opaque raw blocks).
fn write_fast_path(
    out_path: &str,
    header: &Header,
    all: &mut Vec<(i64, Vec<u8>)>,
    workers: NonZeroUsize,
) -> io::Result<()> {
    let t1 = std::time::Instant::now();
    // Stable sort: preserves input order among records sharing the same (refID,pos),
    // matching samtools sort (also stable) so the output is byte-identical.
    all.par_sort_by_key(|(k, _)| *k);
    eprintln!(
        "[bam-insort] fast-path sorted {} reads in RAM in {:.1}s | RSS {:.1} GB",
        all.len(),
        t1.elapsed().as_secs_f64(),
        rss_gb()
    );

    let t2 = std::time::Instant::now();
    let file = File::create(out_path)?;
    let encoder = bgzf::MultithreadedWriter::with_worker_count(workers, file);
    let mut bw = bam::io::Writer::from(encoder);
    bw.write_header(header)?;
    let mut encoder = bw.into_inner();
    for (_, d) in all.iter() {
        write_raw_record(&mut encoder, d)?;
    }
    encoder.finish()?; // flush worker blocks + EOF
    eprintln!(
        "[bam-insort] fast-path wrote {} in {:.1}s | peak RSS {:.1} GB",
        out_path,
        t2.elapsed().as_secs_f64(),
        peak_rss_gb()
    );
    Ok(())
}

/// Bucket path: sort + compress each genomic bucket in parallel, then byte-
/// concatenate the bucket BGZF blobs in genomic order with one EOF marker.
fn write_bucket_path(
    out_path: &str,
    header: &Header,
    buckets: Vec<Bucket>,
    tmp_dir: &str,
    workers: NonZeroUsize,
) -> io::Result<()> {
    let _ = workers;
    // Header as its own BGZF blob (no EOF) -> becomes the start of the file.
    let mut hw = bam::io::Writer::from(bgzf::Writer::new(Vec::new()));
    hw.write_header(header)?;
    let mut hbw = hw.into_inner();
    hbw.flush()?;
    let header_blob = hbw.into_inner();

    // Phase 2: parallel per-bucket sort + compress -> temp BGZF blob (no EOF).
    use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
    let t = std::time::Instant::now();
    let tmp = tmp_dir.to_string();
    let n_buckets = buckets.len();
    let done = AtomicUsize::new(0);
    let recs_done = AtomicU64::new(0);
    let last_decile = AtomicUsize::new(0);
    eprintln!("[bam-insort] phase2 START: sort+compress {n_buckets} genomic buckets in parallel");
    let results: io::Result<Vec<(usize, Option<String>, u64)>> = buckets
        .into_par_iter()
        .enumerate()
        .map(|(idx, b)| {
            let (i, path, nrec) = compress_bucket(idx, b, &tmp)?;
            let c = done.fetch_add(1, Ordering::Relaxed) + 1;
            let r = recs_done.fetch_add(nrec, Ordering::Relaxed) + nrec;
            let decile = c * 10 / n_buckets.max(1);
            if decile > last_decile.load(Ordering::Relaxed)
                && last_decile.fetch_max(decile, Ordering::Relaxed) < decile
            {
                eprintln!(
                    "[bam-insort] phase2 {}0% ({c}/{n_buckets} buckets, {}M reads) | {:.0}s | RSS {:.1} GB",
                    decile,
                    r / 1_000_000,
                    t.elapsed().as_secs_f64(),
                    rss_gb()
                );
            }
            Ok((i, path, nrec))
        })
        .collect();
    let mut blobs: Vec<(usize, String)> = results?
        .into_iter()
        .filter_map(|(i, o, _)| o.map(|p| (i, p)))
        .collect();
    blobs.sort_by_key(|(i, _)| *i);
    eprintln!(
        "[bam-insort] phase2 DONE: {} non-empty buckets in {:.1}s | peak RSS {:.1} GB",
        blobs.len(),
        t.elapsed().as_secs_f64(),
        peak_rss_gb()
    );

    // Phase 3: concatenate header + bucket blobs (genomic order) + single EOF.
    let t = std::time::Instant::now();
    let mut out = BufWriter::new(File::create(out_path)?);
    out.write_all(&header_blob)?;
    for (_, p) in &blobs {
        let mut f = BufReader::new(File::open(p)?);
        io::copy(&mut f, &mut out)?;
    }
    out.write_all(BGZF_EOF)?;
    out.flush()?;
    for (_, p) in &blobs {
        let _ = std::fs::remove_file(p);
    }
    eprintln!("[bam-insort] concatenated {} -> {} in {:.1}s", blobs.len(), out_path, t.elapsed().as_secs_f64());
    Ok(())
}

/// Load one bucket (in-RAM tail + spill file), sort by key, and BGZF-compress to
/// a temp blob WITHOUT an EOF marker (so blobs can be concatenated).
fn compress_bucket(idx: usize, mut b: Bucket, tmp: &str) -> io::Result<(usize, Option<String>, u64)> {
    // Load in input order: spilled records (written in read order across spill
    // rounds) FIRST, then the in-RAM tail (the most recent records). A stable sort
    // then preserves input order among equal (refID,pos), matching samtools sort.
    let mut recs: Vec<(i64, Vec<u8>)> = Vec::new();
    if let Some(sp) = b.spill_path.take() {
        // zstd decoder transparently decodes the concatenated per-round frames.
        let mut r = zstd::stream::read::Decoder::new(BufReader::new(File::open(&sp)?))?;
        let mut d = Vec::new();
        while read_raw_record(&mut r, &mut d)? {
            let key = key_of(&d);
            recs.push((key, std::mem::take(&mut d)));
        }
        let _ = std::fs::remove_file(&sp);
    }
    recs.append(&mut b.mem);
    if recs.is_empty() {
        return Ok((idx, None, 0));
    }
    let nrec = recs.len() as u64;
    recs.sort_by_key(|(k, _)| *k);

    let path = format!("{}/bam-insort.blob.{}.bgz", tmp, idx);
    let mut w = bgzf::Writer::new(BufWriter::new(File::create(&path)?));
    for (_, d) in &recs {
        write_raw_record(&mut w, d)?;
    }
    w.flush()?; // flush final block, NO EOF
    let bufw = w.into_inner();
    bufw.into_inner().map_err(|e| e.into_error())?; // flush BufWriter -> File
    Ok((idx, Some(path), nrec))
}

/// Stamp the @HD line with `SO:coordinate`. If the input had no @HD line, add one.
fn set_coordinate_sorted(header: &mut Header) {
    use noodles_sam::header::record::value::map::{self, header::tag, Map};
    let hd = header
        .header_mut()
        .get_or_insert_with(Map::<map::Header>::default);
    hd.other_fields_mut()
        .insert(tag::SORT_ORDER, bstr::BString::from("coordinate"));
    // Drop any carried-over group-order (e.g. GO:query from the query-grouped input);
    // a coordinate-sorted BAM has no meaningful group order, matching `samtools sort`.
    hd.other_fields_mut().shift_remove(&tag::GROUP_ORDER);
}

/// Build a BAI for the just-written coordinate-sorted BAM (canonical noodles
/// pattern: track BGZF virtual position before/after each record to form the
/// (start, end) chunk and feed the binning `Indexer`).
fn write_bai(bam_path: &str, workers: NonZeroUsize) -> io::Result<()> {
    let file = BufReader::new(File::open(bam_path)?);
    let decoder = bgzf::MultithreadedReader::with_worker_count(workers, file);
    let mut reader = bam::io::Reader::from(decoder);
    let header = reader.read_header()?;

    let mut indexer = Indexer::default();
    let mut record = bam::Record::default();
    let mut start_position = reader.get_ref().virtual_position();

    while reader.read_record(&mut record)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let chunk = Chunk::new(start_position, end_position);

        let alignment_context = match (
            record.reference_sequence_id().transpose()?,
            record.alignment_start().transpose()?,
            record.alignment_end().transpose()?,
        ) {
            (Some(id), Some(start), Some(end)) => {
                let is_mapped = !record.flags().is_unmapped();
                Some((id, start, end, is_mapped))
            }
            _ => None,
        };

        indexer.add_record(alignment_context, chunk)?;
        start_position = end_position;
    }

    let index = indexer.build(header.reference_sequences().len());
    let bai_path = format!("{bam_path}.bai");
    let mut writer = bai::Writer::new(File::create(&bai_path)?);
    writer.write_index(&index)?;
    Ok(())
}

fn num_cpus_env() -> usize {
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(4)
}

// ---- lightweight self-monitoring (Linux /proc) for progress logging ----

/// Current resident set size (GB) from /proc/self/status VmRSS.
fn rss_gb() -> f64 {
    proc_status_kb("VmRSS:") as f64 / (1024.0 * 1024.0)
}
/// Peak resident set size (GB) from /proc/self/status VmHWM.
fn peak_rss_gb() -> f64 {
    proc_status_kb("VmHWM:") as f64 / (1024.0 * 1024.0)
}
fn proc_status_kb(field: &str) -> u64 {
    if let Ok(s) = std::fs::read_to_string("/proc/self/status") {
        for line in s.lines() {
            if let Some(rest) = line.strip_prefix(field) {
                if let Some(kb) = rest.split_whitespace().next().and_then(|v| v.parse::<u64>().ok()) {
                    return kb;
                }
            }
        }
    }
    0
}

/// Total CPU seconds (utime+stime) used by this process, from /proc/self/stat.
fn cpu_secs() -> f64 {
    if let Ok(s) = std::fs::read_to_string("/proc/self/stat") {
        // fields after the "(comm)" may contain spaces; split on the last ')'.
        if let Some(idx) = s.rfind(')') {
            let rest: Vec<&str> = s[idx + 1..].split_whitespace().collect();
            // rest[0]=state (field 3); utime=field14 -> rest[11], stime=field15 -> rest[12]
            if rest.len() > 12 {
                let ut = rest[11].parse::<f64>().unwrap_or(0.0);
                let st = rest[12].parse::<f64>().unwrap_or(0.0);
                return (ut + st) / 100.0; // USER_HZ = 100 on Linux
            }
        }
    }
    0.0
}

/// Tracks avg CPU cores used between successive calls (delta CPU / delta wall).
struct ProcMon {
    last_wall: std::time::Instant,
    last_cpu: f64,
}
impl ProcMon {
    fn new() -> Self {
        ProcMon {
            last_wall: std::time::Instant::now(),
            last_cpu: cpu_secs(),
        }
    }
    fn cores(&mut self) -> f64 {
        let w = std::time::Instant::now();
        let c = cpu_secs();
        let dt = (w - self.last_wall).as_secs_f64().max(1e-6);
        let cores = (c - self.last_cpu) / dt;
        self.last_wall = w;
        self.last_cpu = c;
        cores
    }
}

/// Read MemAvailable (kB) from /proc/meminfo; fall back to a conservative 8 GB.
fn mem_available_bytes() -> u64 {
    if let Ok(s) = std::fs::read_to_string("/proc/meminfo") {
        for line in s.lines() {
            if let Some(rest) = line.strip_prefix("MemAvailable:") {
                if let Some(kb) = rest.split_whitespace().next().and_then(|v| v.parse::<u64>().ok()) {
                    return kb * 1024;
                }
            }
        }
    }
    8 * 1024 * 1024 * 1024
}
