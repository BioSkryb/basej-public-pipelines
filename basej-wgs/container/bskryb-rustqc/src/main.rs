// bskryb-qc: single-pass DNA QC metrics from a coordinate-sorted BAM.
//
// Phase 1a rewrite of the Sentieon/Picard AlignmentStat + InsertSizeMetricAlgo
// outputs consumed by the basej-dnaqc metrics Parquet. Follows the rewrites.bio
// principle "emulate exactly": we reproduce the interval-restricted per-category
// counts Picard/Sentieon produce, and validate field-by-field against the original.
//
// Credits (tools being emulated):
//   - Sentieon driver AlignmentStat / InsertSizeMetricAlgo (Sentieon Inc.)
//   - Picard CollectAlignmentSummaryMetrics / CollectInsertSizeMetrics (Broad Institute)
//
// AI provenance: implementation is AI-assisted; correctness is determined by output
// comparison against the originals on real BioSkryb biosamples, not code review alone.

use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

use noodles_bam as bam;
use noodles_sam::alignment::record::data::field::Tag;

mod preseq;

// ---------------------------------------------------------------------------
// Intervals: contig name -> sorted list of (start0, end0) half-open, 0-based.
// ---------------------------------------------------------------------------
#[derive(Clone)]
struct Intervals {
    /// Intervals as read (sorted). Hashed for the GC-window index and used by the
    /// coverage-subset check, so keep it exactly as before.
    map: HashMap<String, Vec<(i64, i64)>>,
    /// Same union, merged into disjoint sorted runs: binary-search lookups only.
    merged: HashMap<String, Vec<(i64, i64)>>,
    enabled: bool,
    /// --interval-overlap: admit a read when its aligned span OVERLAPS an interval
    /// (Sentieon/Picard `--interval` semantics), instead of requiring its start position to
    /// lie inside one. Matters for small targets (exome); WGS intervals are chromosome-scale
    /// so the two rules agree there and the default stays start-in-interval.
    overlap: bool,
}

/// Report the first coverage interval that is not inside the base-metrics set, if any.
///
/// Both hot loops gate on the base set first and only then consult the coverage set, so a
/// coverage region outside the base set would be silently excluded from coverage rather than
/// counted. Sentieon's own split has coverage as the narrower set, so this is a config error
/// rather than a case to support -- fail loudly instead of reporting quietly-wrong coverage.
fn coverage_subset_violation(cov: &Intervals, base: &Intervals) -> Option<String> {
    if !base.enabled {
        return None; // no base restriction, so everything is already in range
    }
    if !cov.enabled {
        return Some("coverage interval set is unrestricted while --intervals restricts".into());
    }
    for (contig, ivs) in &cov.map {
        let base_ivs = match base.map.get(contig) {
            Some(v) => v,
            None => {
                return Some(format!(
                    "contig {contig} is in --coverage-intervals but absent from --intervals"
                ))
            }
        };
        for &(s, e) in ivs {
            if !base_ivs.iter().any(|&(bs, be)| s >= bs && e <= be) {
                return Some(format!(
                    "coverage interval {contig}:{}-{} is not contained in any --intervals region",
                    s + 1, e
                ));
            }
        }
    }
    None
}

impl Intervals {
    fn none() -> Self {
        Intervals { map: HashMap::new(), merged: HashMap::new(), enabled: false, overlap: false }
    }

    fn from_bed(path: &str) -> std::io::Result<Self> {
        // Supports both BED (0-based, half-open) and Picard/GATK interval_list
        // (1-based, inclusive, with @-prefixed SAM headers). The interval_list
        // form is what genomes.config's `wgs_or_target_intervals` uses, so the
        // WGS coverage regions parse identically to what Sentieon consumed.
        let is_interval_list = path.ends_with(".interval_list") || path.ends_with(".intervals")
            || path.ends_with(".interval_list.gz");
        let mut map: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
        let f = BufReader::new(File::open(path)?);
        for line in f.lines() {
            let line = line?;
            if line.is_empty() || line.starts_with('#') || line.starts_with("track")
                || line.starts_with('@') {
                continue;
            }
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 3 {
                continue;
            }
            // interval_list is 1-based inclusive -> convert to 0-based half-open.
            let (start, end): (i64, i64) = if is_interval_list {
                (cols[1].parse::<i64>().unwrap_or(1) - 1, cols[2].parse().unwrap_or(0))
            } else {
                (cols[1].parse().unwrap_or(0), cols[2].parse().unwrap_or(0))
            };
            map.entry(cols[0].to_string()).or_default().push((start, end));
        }
        for v in map.values_mut() {
            v.sort_unstable();
        }
        // Lookup copy: overlapping/abutting intervals merged. The union is unchanged, so
        // point membership is identical to the old linear scan; disjoint sorted runs allow
        // binary search (exome panels have ~200k targets).
        let mut merged: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
        for (c, v) in &map {
            let mut m: Vec<(i64, i64)> = Vec::with_capacity(v.len());
            for &(s, e) in v {
                match m.last_mut() {
                    Some(last) if s <= last.1 => { if e > last.1 { last.1 = e; } }
                    _ => m.push((s, e)),
                }
            }
            merged.insert(c.clone(), m);
        }
        Ok(Intervals { map, merged, enabled: true, overlap: false })
    }

    // pos is 1-based alignment start; test overlap of the single base at pos.
    fn contains(&self, contig: &str, pos1: i64) -> bool {
        if !self.enabled {
            return true;
        }
        if pos1 <= 0 {
            return false;
        }
        let p0 = pos1 - 1; // to 0-based
        match self.merged.get(contig) {
            None => false,
            Some(ivs) => {
                // last interval starting at or before p0 (list is disjoint + sorted)
                let i = ivs.partition_point(|&(s, _)| s <= p0);
                i > 0 && p0 < ivs[i - 1].1
            }
        }
    }

    /// Placement test used by the alignment/GC/insert accumulators. Start-in-interval by
    /// default; with `overlap`, any overlap of the 1-based inclusive span [start1, end1].
    fn admits(&self, contig: &str, start1: i64, end1: i64) -> bool {
        if !self.overlap {
            return self.contains(contig, start1);
        }
        if !self.enabled {
            return true;
        }
        if start1 <= 0 {
            return false;
        }
        let (a0, b0) = (start1 - 1, end1.max(start1)); // 0-based half-open [a0, b0)
        match self.merged.get(contig) {
            None => false,
            Some(ivs) => {
                // first interval ending after a0; it overlaps iff it starts before b0
                let i = ivs.partition_point(|&(_, e)| e <= a0);
                i < ivs.len() && ivs[i].0 < b0
            }
        }
    }
}

/// 1-based inclusive alignment end from the CIGAR (M/=/X/D/N consume reference).
fn alignment_end1(record: &bam::Record, start1: i64) -> i64 {
    use noodles_sam::alignment::record::cigar::op::Kind;
    let mut span: i64 = 0;
    for op in record.cigar().iter().flatten() {
        match op.kind() {
            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch
            | Kind::Deletion | Kind::Skip => span += op.len() as i64,
            _ => {}
        }
    }
    start1 + span.max(1) - 1
}

// ---------------------------------------------------------------------------
// Chimera feature cube (modular substrate for chimera composition).
// Axes: LOCUS(2) x ORIENT(4) x DIST(6) x SA(2) = 96 cells. Every aligned-in-pair
// read is binned into exactly one cell. Named chimera mechanism classes are
// declarative sums over cube cells (see chimera_classes.yaml), so new/experimental
// classes can be defined WITHOUT re-running alignment. The set of "non-proper" cells
// reproduces the union PCT_CHIMERAS exactly (proper = same_contig & proper_FR &
// dist<100kb & no_SA). Discovery: cluster the cube across cells/chemistries to find
// recurrent structures not yet named, then add a rule.
// ---------------------------------------------------------------------------
const CUBE_LOCUS: usize = 2; //  0=same_contig 1=diff_contig
const CUBE_ORIENT: usize = 4; // 0=proper_FR 1=FF 2=RR 3=everted
const CUBE_DIST: usize = 6; //   0=na(diff) 1=<1kb 2=1-10kb 3=10-100kb 4=100kb-1Mb 5=>1Mb
const CUBE_SA: usize = 2; //     0=no_SA 1=has_SA
const N_CHIM_CELLS: usize = CUBE_LOCUS * CUBE_ORIENT * CUBE_DIST * CUBE_SA; // 96

#[derive(Clone)]
struct ChimCube([u64; N_CHIM_CELLS]);
impl Default for ChimCube {
    fn default() -> Self { ChimCube([0u64; N_CHIM_CELLS]) }
}

#[inline]
fn cube_idx(locus: usize, orient: usize, dist: usize, sa: usize) -> usize {
    ((locus * CUBE_ORIENT + orient) * CUBE_DIST + dist) * CUBE_SA + sa
}

/// Bin one aligned-in-pair read into a chimera-cube cell index.
#[inline]
fn chim_cube_bin(same_contig: bool, fr_proper: bool, read_rev: bool, mate_rev: bool,
                 tlen: i64, has_sa: bool) -> usize {
    let locus = if same_contig { 0 } else { 1 };
    // proper_FR is only meaningful within a contig; inter-chromosomal reads have TLEN=0
    // (which would spuriously satisfy fr_proper), so classify them by strand combo.
    let orient = if same_contig && fr_proper { 0 }
        else if !read_rev && !mate_rev { 1 } // FF (same-strand)
        else if read_rev && mate_rev { 2 }   // RR (same-strand)
        else { 3 };                          // everted / opposite-strand (not proper_FR)
    let dist = if !same_contig { 0 } else {
        let a = tlen.abs();
        if a < 1_000 { 1 } else if a < 10_000 { 2 } else if a < 100_000 { 3 }
        else if a < 1_000_000 { 4 } else { 5 }
    };
    let sa = usize::from(has_sa);
    cube_idx(locus, orient, dist, sa)
}

// ---------------------------------------------------------------------------
// Per-category accumulators (FIRST_OF_PAIR, SECOND_OF_PAIR).
// ---------------------------------------------------------------------------
#[derive(Default, Clone)]
struct Cat {
    total_reads: u64,
    pf_reads_aligned: u64,
    read_length_sum: u64,
    reads_aligned_in_pairs: u64,
    neg_strand_aligned: u64,
    mismatch_sum: u64,   // MD-based substitutions (exact, excludes indels)
    aligned_bases: u64,  // CIGAR M/=/X
    indel_events: u64,   // number of I/D CIGAR operations
    indel_bases: u64,    // I/D base lengths
    // chimera total + component breakdown (a read may satisfy several; total = union)
    chimeric_reads: u64,
    chim_diff_contig: u64,
    chim_bad_orient: u64,
    chim_large_insert: u64,
    chim_split_sa: u64,
    // bad_orientation sub-types (mate on same contig, non-FR): FF/RR = same-strand
    // (foldback / inverted-repeat WGA junctions), RF = everted "outie" (tandem-dup /
    // circular / read-through). FF+RR+RF == chim_bad_orient.
    chim_bo_ff: u64,
    chim_bo_rr: u64,
    chim_bo_rf: u64,
    // bad_orientation intra-pair distance buckets (|TLEN|): <1kb dominated => local
    // hairpin junctions; larger => distal rearrangements. Sum == chim_bad_orient.
    chim_bo_dist_lt1kb: u64,
    chim_bo_dist_1_10kb: u64,
    chim_bo_dist_10_100kb: u64,
    chim_bo_dist_gt100kb: u64,
    // modular chimera feature cube (LOCUS x ORIENT x DIST x SA)
    chim_cube: ChimCube,
    // high-quality (MAPQ >= 20) sub-metrics
    hq_reads: u64,
    hq_aligned_bases: u64,
    hq_q20_bases: u64,
    hq_mismatch_sum: u64,
    // Single-end (UNPAIRED) chimera denominator: mapped reads with MAPQ >= 20, as Picard
    // CollectAlignmentSummaryMetrics counts for fragments. Always 0 for paired categories,
    // whose denominator stays reads_aligned_in_pairs (paired output unchanged).
    chim_fragment_den: u64,
}

impl Cat {
    fn add(&mut self, o: &Cat) {
        self.total_reads += o.total_reads;
        self.pf_reads_aligned += o.pf_reads_aligned;
        self.read_length_sum += o.read_length_sum;
        self.reads_aligned_in_pairs += o.reads_aligned_in_pairs;
        self.neg_strand_aligned += o.neg_strand_aligned;
        self.mismatch_sum += o.mismatch_sum;
        self.aligned_bases += o.aligned_bases;
        self.indel_events += o.indel_events;
        self.indel_bases += o.indel_bases;
        self.chimeric_reads += o.chimeric_reads;
        self.chim_diff_contig += o.chim_diff_contig;
        self.chim_bad_orient += o.chim_bad_orient;
        self.chim_large_insert += o.chim_large_insert;
        self.chim_split_sa += o.chim_split_sa;
        self.chim_bo_ff += o.chim_bo_ff;
        self.chim_bo_rr += o.chim_bo_rr;
        self.chim_bo_rf += o.chim_bo_rf;
        self.chim_bo_dist_lt1kb += o.chim_bo_dist_lt1kb;
        self.chim_bo_dist_1_10kb += o.chim_bo_dist_1_10kb;
        self.chim_bo_dist_10_100kb += o.chim_bo_dist_10_100kb;
        self.chim_bo_dist_gt100kb += o.chim_bo_dist_gt100kb;
        for i in 0..N_CHIM_CELLS {
            self.chim_cube.0[i] += o.chim_cube.0[i];
        }
        self.hq_reads += o.hq_reads;
        self.hq_aligned_bases += o.hq_aligned_bases;
        self.hq_q20_bases += o.hq_q20_bases;
        self.hq_mismatch_sum += o.hq_mismatch_sum;
        self.chim_fragment_den += o.chim_fragment_den;
    }
}

/// Sub-classify a bad_orientation read (same contig, non-FR): strand combination
/// (FF/RR = same-strand foldback/inverted; RF = everted "outie") + intra-pair
/// distance bucket (|TLEN|). FF+RR+RF and the 4 distance buckets each sum to
/// chim_bad_orient. Called from both the single-thread and parallel accumulation paths.
#[inline]
fn classify_bad_orient(cat: &mut Cat, read_rev: bool, mate_rev: bool, tlen: i64) {
    if !read_rev && !mate_rev {
        cat.chim_bo_ff += 1;
    } else if read_rev && mate_rev {
        cat.chim_bo_rr += 1;
    } else {
        cat.chim_bo_rf += 1;
    }
    let a = tlen.abs();
    if a < 1_000 {
        cat.chim_bo_dist_lt1kb += 1;
    } else if a < 10_000 {
        cat.chim_bo_dist_1_10kb += 1;
    } else if a < 100_000 {
        cat.chim_bo_dist_10_100kb += 1;
    } else {
        cat.chim_bo_dist_gt100kb += 1;
    }
}

/// Count substitution bases encoded in an MD tag (excludes deletions `^...`).
fn md_mismatches(md: &str) -> u64 {
    let b = md.as_bytes();
    let mut i = 0;
    let mut n = 0u64;
    while i < b.len() {
        let c = b[i];
        if c.is_ascii_digit() {
            i += 1;
        } else if c == b'^' {
            i += 1;
            while i < b.len() && b[i].is_ascii_alphabetic() {
                i += 1;
            }
        } else if c.is_ascii_alphabetic() {
            n += 1;
            i += 1;
        } else {
            i += 1;
        }
    }
    n
}

/// Reference-consuming length of a CIGAR string (from the MC mate-cigar tag).
fn ref_len_from_cigar_str(cigar: &str) -> i64 {
    let mut len = 0i64;
    let mut num = 0i64;
    for c in cigar.bytes() {
        if c.is_ascii_digit() {
            num = num * 10 + (c - b'0') as i64;
        } else {
            if matches!(c, b'M' | b'D' | b'N' | b'=' | b'X') {
                len += num;
            }
            num = 0;
        }
    }
    len
}

// ---------------------------------------------------------------------------
// GC bias (Picard CollectGcBiasMetrics-equivalent). Requires an indexed FASTA.
// ---------------------------------------------------------------------------
struct FaiEntry {
    len: u64,
    offset: u64,
    linebases: u64,
    linewidth: u64,
}

fn load_fai(path: &str) -> std::io::Result<HashMap<String, FaiEntry>> {
    let mut map = HashMap::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() >= 5 {
            map.insert(
                f[0].to_string(),
                FaiEntry {
                    len: f[1].parse().unwrap_or(0),
                    offset: f[2].parse().unwrap_or(0),
                    linebases: f[3].parse().unwrap_or(0),
                    linewidth: f[4].parse().unwrap_or(0),
                },
            );
        }
    }
    Ok(map)
}

/// Load a contig's bases (uppercased, newlines stripped) using its .fai entry.
fn load_contig(fa: &mut File, e: &FaiEntry) -> std::io::Result<Vec<u8>> {
    use std::io::{Read, Seek, SeekFrom};
    let n_lines = if e.linebases > 0 { (e.len + e.linebases - 1) / e.linebases } else { 0 };
    let raw_len = if e.len == 0 { 0 } else { (n_lines - 1) * e.linewidth + (e.len - (n_lines - 1) * e.linebases) };
    fa.seek(SeekFrom::Start(e.offset))?;
    let mut raw = vec![0u8; raw_len as usize];
    fa.read_exact(&mut raw)?;
    let mut bases = Vec::with_capacity(e.len as usize);
    for b in raw {
        if b != b'\n' && b != b'\r' {
            bases.push(b.to_ascii_uppercase());
        }
    }
    Ok(bases)
}

/// GC value (0..window, via integer division *100/window) of the window [start,start+window)
/// over 0-based bases, or -1 if the window has > 4 Ns or runs off the end. Matches
/// Picard GcBiasUtils.calculateGc semantics for a single window.
fn window_gc_at(bases: &[u8], start: usize, window: usize) -> i32 {
    if start + window > bases.len() {
        return 0; // Picard default for uncomputed indices
    }
    let mut gc_count: i32 = 0;
    let mut n_count: i32 = 0;
    for &b in &bases[start..start + window] {
        if b == b'G' || b == b'C' {
            gc_count += 1;
        } else if b == b'N' {
            n_count += 1;
        }
    }
    if n_count > 4 {
        -1
    } else {
        (gc_count * 100) / window as i32
    }
}

/// Accumulate windowsByGc over a contig's bases (streaming; no full gc[] array).
/// Matches Picard calculateRefWindowsByGc: windows start at i in 1..(len-window),
/// integer-division GC, skip windows with > 4 Ns.
fn add_windows_by_gc(bases: &[u8], window: usize, hist: &mut [u64; 101]) {
    let n = bases.len();
    if n <= window {
        return;
    }
    let last_window_start = n - window;
    let is_gc = |b: u8| b == b'G' || b == b'C';
    let mut gc_count: i32 = 0;
    let mut n_count: i32 = 0;
    for i in 1..last_window_start {
        if i == 1 {
            for j in 1..(1 + window) {
                let b = bases[j];
                if is_gc(b) { gc_count += 1; } else if b == b'N' { n_count += 1; }
            }
        } else {
            let newb = bases[i + window - 1];
            if is_gc(newb) { gc_count += 1; } else if newb == b'N' { n_count += 1; }
            let oldb = bases[i - 1];
            if is_gc(oldb) { gc_count -= 1; } else if oldb == b'N' { n_count -= 1; }
        }
        if n_count <= 4 {
            hist[((gc_count * 100) / window as i32) as usize] += 1;
        }
    }
}

/// Stable FNV-1a hash of the normalized interval set (sorted "contig\tstart\tend" lines),
/// used to key/validate a cached GC-window index against the intervals in use.
fn intervals_hash(intervals: &Intervals) -> String {
    let mut lines: Vec<String> = Vec::new();
    for (c, ivs) in &intervals.map {
        for (s, e) in ivs {
            lines.push(format!("{c}\t{s}\t{e}"));
        }
    }
    lines.sort();
    let mut h: u64 = 0xcbf29ce484222325;
    for l in lines {
        for b in l.bytes() {
            h ^= b as u64;
            h = h.wrapping_mul(0x100000001b3);
        }
        h ^= b'\n' as u64;
        h = h.wrapping_mul(0x100000001b3);
    }
    format!("{h:016x}")
}

/// Write a GC-window index file (windowsByGc + metadata) for caching.
fn write_gc_index(path: &str, window: usize, ivhash: &str, hist: &[u64; 101]) -> std::io::Result<()> {
    let mut f = File::create(path)?;
    writeln!(f, "# bskryb-qc gc-index v1")?;
    writeln!(f, "# window\t{window}")?;
    writeln!(f, "# intervals_sha\t{ivhash}")?;
    writeln!(f, "GC\tWINDOWS")?;
    for (gc, &c) in hist.iter().enumerate() {
        writeln!(f, "{gc}\t{c}")?;
    }
    Ok(())
}

/// Read a cached GC-window index; returns (windowsByGc, window, intervals_sha).
fn read_gc_index(path: &str) -> std::io::Result<([u64; 101], usize, String)> {
    let mut hist = [0u64; 101];
    let mut window = 0usize;
    let mut ivhash = String::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix("# window\t") {
            window = rest.trim().parse().unwrap_or(0);
        } else if let Some(rest) = line.strip_prefix("# intervals_sha\t") {
            ivhash = rest.trim().to_string();
        } else if line.starts_with('#') || line.starts_with("GC\t") {
            continue;
        } else {
            let f: Vec<&str> = line.split('\t').collect();
            if f.len() == 2 {
                if let (Ok(gc), Ok(c)) = (f[0].parse::<usize>(), f[1].parse::<u64>()) {
                    if gc < 101 { hist[gc] = c; }
                }
            }
        }
    }
    Ok((hist, window, ivhash))
}

struct GcBias {
    window: usize,
    fa_path: String,
    fai: HashMap<String, FaiEntry>,
    windows_by_gc: [u64; 101],
    reads_by_gc: [u64; 101],
    total_clusters: u64,
    total_aligned_reads: u64,
    cur_contig: String,
    bases: Vec<u8>,          // current contig bases (loaded on demand)
    last_window_start: usize, // for the current contig
}

impl GcBias {
    /// `windows_cache`: if Some, use the precomputed windowsByGc (single-pass: no
    /// full-genome scan). If None, compute it over `contigs` (standalone fallback).
    fn new(
        fa_path: &str,
        window: usize,
        contigs: &[String],
        windows_cache: Option<[u64; 101]>,
    ) -> std::io::Result<Self> {
        let fai = load_fai(&format!("{fa_path}.fai"))?;
        let mut gb = GcBias {
            window,
            fa_path: fa_path.to_string(),
            fai,
            windows_by_gc: [0u64; 101],
            reads_by_gc: [0u64; 101],
            total_clusters: 0,
            total_aligned_reads: 0,
            cur_contig: String::new(),
            bases: Vec::new(),
            last_window_start: 0,
        };
        match windows_cache {
            Some(w) => gb.windows_by_gc = w,
            None => {
                // Fallback: compute windowsByGc over the interval contigs.
                let mut fa = File::open(fa_path)?;
                for c in contigs {
                    if let Some(e) = gb.fai.get(c) {
                        let bases = load_contig(&mut fa, e)?;
                        add_windows_by_gc(&bases, window, &mut gb.windows_by_gc);
                    }
                }
            }
        }
        Ok(gb)
    }

    fn ensure_contig(&mut self, contig: &str) {
        if self.cur_contig == contig {
            return;
        }
        self.bases = Vec::new();
        self.last_window_start = 0;
        if let Some(e) = self.fai.get(contig) {
            if let Ok(mut fa) = File::open(&self.fa_path) {
                if let Ok(bases) = load_contig(&mut fa, e) {
                    self.last_window_start = if bases.len() > self.window { bases.len() - self.window } else { 0 };
                    self.bases = bases;
                }
            }
        }
        self.cur_contig = contig.to_string();
    }

    /// Add a mapped read. `start1`/`end1` are 1-based alignment start/end (inclusive).
    /// Per-read window GC is computed on the fly (Picard indexes gc[pos] with the
    /// 1-based coordinate used directly as a 0-based offset — replicated here).
    fn add_read(&mut self, contig: &str, start1: i64, end1: i64, reverse: bool, first_of_pair: bool) {
        self.ensure_contig(contig);
        if first_of_pair {
            self.total_clusters += 1;
        }
        self.total_aligned_reads += 1;
        let pos = if reverse { end1 - self.window as i64 } else { start1 };
        if pos <= 0 || self.bases.is_empty() {
            return;
        }
        let p = pos as usize;
        // Picard: gc array only populated for 1..lastWindowStart; elsewhere default 0.
        let windowgc = if p >= 1 && p < self.last_window_start {
            window_gc_at(&self.bases, p, self.window)
        } else {
            0
        };
        if windowgc >= 0 {
            self.reads_by_gc[windowgc as usize] += 1;
        }
    }

    /// Merge another GcBias (same reference/windows) into this one. reads_by_gc and
    /// the read/cluster counts are additive across region-sharded workers;
    /// windows_by_gc is identical (from the shared cache) so it is left as-is.
    fn merge(&mut self, o: &GcBias) {
        for i in 0..101 {
            self.reads_by_gc[i] += o.reads_by_gc[i];
        }
        self.total_clusters += o.total_clusters;
        self.total_aligned_reads += o.total_aligned_reads;
    }

    fn norm_cov(&self, mean_rpw: f64, start: usize, end: usize) -> f64 {
        let mut windows_total: u64 = 0;
        let mut sum: f64 = 0.0;
        for i in start..=end {
            if self.windows_by_gc[i] != 0 {
                sum += self.reads_by_gc[i] as f64;
                windows_total += self.windows_by_gc[i];
            }
        }
        if windows_total == 0 { 0.0 } else { sum / (windows_total as f64 * mean_rpw) }
    }
}



/// preseq gc_extrap coverage-count histogram builder: tiles the genome into `bin_size`
/// (10bp) bins and counts reads per bin via a coordinate sweep, then histograms the
/// per-bin coverage counts. Deterministic binning: a read is registered on a bin if it
/// covers >= bin_size/2 bases of it (round-at-0.5, matching preseq's expected binning).
struct PreseqAcc {
    bin_size: i64,
    hist: Vec<f64>,
    cur_ref: i64,
    pending: std::collections::BTreeMap<i64, u32>,
    rng: rand_chacha::ChaCha8Rng,
}

impl PreseqAcc {
    fn new(bin_size: i64) -> Self {
        use rand::SeedableRng;
        PreseqAcc {
            bin_size,
            hist: vec![0.0; 2],
            cur_ref: -1,
            pending: std::collections::BTreeMap::new(),
            rng: rand_chacha::ChaCha8Rng::seed_from_u64(408),
        }
    }

    fn finalize_bin(&mut self, count: u32) {
        let c = count as usize;
        if self.hist.len() < c + 1 {
            self.hist.resize(c + 1, 0.0);
        }
        self.hist[c] += 1.0;
    }

    fn flush_before(&mut self, bin_idx: i64) {
        let keys: Vec<i64> = self.pending.range(..bin_idx).map(|(k, _)| *k).collect();
        for k in keys {
            if let Some(c) = self.pending.remove(&k) {
                self.finalize_bin(c);
            }
        }
    }

    fn flush_all(&mut self) {
        let vals: Vec<u32> = self.pending.values().copied().collect();
        for c in vals {
            self.finalize_bin(c);
        }
        self.pending.clear();
    }

    fn add_read(&mut self, ref_id: i64, start1: i64, ref_span: i64, ov_s: i64, ov_e: i64) {
        use rand::Rng;
        if ref_id < 0 || ref_span <= 0 {
            return;
        }
        let start0 = start1 - 1;
        let end0 = start0 + ref_span; // half-open
        let first_bin = start0 / self.bin_size;
        if ref_id != self.cur_ref {
            self.flush_all();
            self.cur_ref = ref_id;
        } else {
            self.flush_before(first_bin);
        }
        // Probabilistic per-bin emission (matches preseq SplitMappedRead). Bases already
        // covered by this read's mate (overlap region [ov_s,ov_e)) are excluded so each
        // fragment's overlap counts once (matches bam2mr fragment merging).
        let last_bin = (end0 - 1) / self.bin_size;
        for b in first_bin..=last_bin {
            let bin_start = b * self.bin_size;
            let bin_end = bin_start + self.bin_size;
            let lo = start0.max(bin_start);
            let hi = end0.min(bin_end);
            let mut covered = hi - lo;
            // subtract the part of this bin already covered by the mate
            if ov_e > ov_s {
                let olo = lo.max(ov_s);
                let ohi = hi.min(ov_e);
                if ohi > olo {
                    covered -= ohi - olo;
                }
            }
            if covered <= 0 {
                continue;
            }
            let frac = covered as f64 / self.bin_size as f64;
            if self.rng.gen::<f64>() <= frac {
                *self.pending.entry(b).or_insert(0) += 1;
            }
        }
    }
}


/// coverage cap) with a per-filter excluded-base breakdown. Requires a coordinate-sorted
/// BAM. GENOME_TERRITORY (non-N reference bases) is supplied by the caller.
struct CovAcc {
    enabled: bool,
    territory: u64,
    cap: u32,
    hist: Vec<u64>,
    covered_positions: u64,
    counted_bases: u64,
    exc_capped: u64,
    cur_ref: i64,
    pending: std::collections::BTreeMap<i64, u32>,
    exc_mapq: u64,
    exc_dupe: u64,
    exc_unpaired: u64,
    exc_baseq: u64,
    exc_overlap: u64,
}

impl CovAcc {
    fn new(territory: u64, cap: u32, enabled: bool) -> Self {
        CovAcc {
            enabled,
            territory,
            cap,
            hist: vec![0u64; cap as usize + 1],
            covered_positions: 0,
            counted_bases: 0,
            exc_capped: 0,
            cur_ref: -1,
            pending: std::collections::BTreeMap::new(),
            exc_mapq: 0,
            exc_dupe: 0,
            exc_unpaired: 0,
            exc_baseq: 0,
            exc_overlap: 0,
        }
    }

    fn finalize_pos(&mut self, depth: u32) {
        let capped = depth.min(self.cap);
        if depth > self.cap {
            self.exc_capped += (depth - self.cap) as u64;
        }
        self.hist[capped as usize] += 1;
        if depth >= 1 {
            self.covered_positions += 1;
        }
    }

    fn flush_before(&mut self, up_to: i64) {
        let keys: Vec<i64> = self.pending.range(..up_to).map(|(k, _)| *k).collect();
        for k in keys {
            if let Some(d) = self.pending.remove(&k) {
                self.finalize_pos(d);
            }
        }
    }

    fn flush_all(&mut self) {
        let vals: Vec<u32> = self.pending.values().copied().collect();
        for d in vals {
            self.finalize_pos(d);
        }
        self.pending.clear();
    }

    fn on_new_read(&mut self, ref_id: i64, start0: i64) {
        if ref_id != self.cur_ref {
            self.flush_all();
            self.cur_ref = ref_id;
        } else {
            self.flush_before(start0);
        }
    }

    fn add_position(&mut self, pos0: i64) {
        *self.pending.entry(pos0).or_insert(0) += 1;
        self.counted_bases += 1;
    }

    fn mean(&self) -> f64 {
        if self.territory == 0 { return 0.0; }
        (self.counted_bases - self.exc_capped) as f64 / self.territory as f64
    }

    fn pct_ge(&self, n: usize) -> f64 {
        if self.territory == 0 { return 0.0; }
        let ge: u64 = self.hist.iter().enumerate().filter(|(d, _)| *d >= n).map(|(_, c)| *c).sum();
        ge as f64 / self.territory as f64
    }

    fn median(&self) -> f64 {
        let zeros = self.territory.saturating_sub(self.covered_positions);
        let half = self.territory / 2;
        let mut cum = zeros;
        if cum > half { return 0.0; }
        for (d, c) in self.hist.iter().enumerate().skip(1) {
            cum += *c;
            if cum > half { return d as f64; }
        }
        0.0
    }

    /// Median absolute deviation of per-position coverage from the median
    /// (Picard WgsMetrics MAD_COVERAGE). Deviations are integer (median is an
    /// integer bin); the deviation distribution is built from the coverage
    /// histogram (including zero-coverage territory) and its median taken with
    /// the same `cum > half` rule as `median()`.
    fn mad(&self, median: f64) -> f64 {
        if self.territory == 0 { return 0.0; }
        let zeros = self.territory.saturating_sub(self.covered_positions);
        let mut dev_counts: std::collections::BTreeMap<u64, u64> = std::collections::BTreeMap::new();
        if zeros > 0 {
            let dev = median.abs().round() as u64; // |0 - median|
            *dev_counts.entry(dev).or_insert(0) += zeros;
        }
        for (d, c) in self.hist.iter().enumerate().skip(1) {
            if *c == 0 { continue; }
            let dev = (d as f64 - median).abs().round() as u64;
            *dev_counts.entry(dev).or_insert(0) += *c;
        }
        let half = self.territory / 2;
        let mut cum = 0u64;
        for (dev, c) in dev_counts.iter() {
            cum += *c;
            if cum > half { return *dev as f64; }
        }
        0.0
    }

    fn sd(&self, mean: f64) -> f64 {
        if self.territory == 0 { return 0.0; }
        let zeros = self.territory.saturating_sub(self.covered_positions);
        let mut ss = mean * mean * zeros as f64;
        for (d, c) in self.hist.iter().enumerate().skip(1) {
            let diff = d as f64 - mean;
            ss += diff * diff * (*c as f64);
        }
        (ss / self.territory as f64).sqrt()
    }

    /// Merge another (fully-flushed) CovAcc into this one. Used by the parallel
    /// region-sharded path — each worker accumulates a disjoint set of contigs and
    /// the histograms/exclusion counters are additive. `pending` must be empty
    /// (workers flush at end of each contig / at finish).
    fn merge(&mut self, o: &CovAcc) {
        debug_assert!(o.pending.is_empty());
        for (i, c) in o.hist.iter().enumerate() {
            self.hist[i] += *c;
        }
        self.covered_positions += o.covered_positions;
        self.counted_bases += o.counted_bases;
        self.exc_capped += o.exc_capped;
        self.exc_mapq += o.exc_mapq;
        self.exc_dupe += o.exc_dupe;
        self.exc_unpaired += o.exc_unpaired;
        self.exc_baseq += o.exc_baseq;
        self.exc_overlap += o.exc_overlap;
    }
}

fn f6(v: f64) -> String {
    format!("{:.6}", v)
}

/// Higher precision for Lorenz/Gini parity comparisons against the original tool.
fn f9(v: f64) -> String {
    format!("{:.12}", v)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let result = if args.len() > 1 && args[1] == "gc-index" {
        gc_index(&args)
    } else {
        run()
    };
    if let Err(e) = result {
        eprintln!("bskryb-qc error: {e}");
        std::process::exit(1);
    }
}

/// Subcommand: precompute the windowsByGc index for a reference + interval set.
/// Writes a small cacheable artifact (used at runtime via --gc-windows).
fn gc_index(args: &[String]) -> std::io::Result<()> {
    let mut reference = String::new();
    let mut bed_path = String::new();
    let mut window: usize = 100;
    let mut out_path = String::new();
    let mut i = 2;
    while i < args.len() {
        match args[i].as_str() {
            "--reference" => { reference = args[i + 1].clone(); i += 2; }
            "--intervals" => { bed_path = args[i + 1].clone(); i += 2; }
            "--window" => { window = args[i + 1].parse().unwrap_or(100); i += 2; }
            "--out" => { out_path = args[i + 1].clone(); i += 2; }
            other => { eprintln!("unknown arg: {other}"); std::process::exit(2); }
        }
    }
    if reference.is_empty() || out_path.is_empty() {
        eprintln!("usage: bskryb-qc gc-index --reference <genome.fa> [--intervals <bed>] [--window 100] --out <index.txt>");
        std::process::exit(2);
    }
    let intervals = if bed_path.is_empty() { Intervals::none() } else { Intervals::from_bed(&bed_path)? };
    let ivhash = if intervals.enabled { intervals_hash(&intervals) } else { "none".to_string() };

    let fai = load_fai(&format!("{reference}.fai"))?;
    let mut contigs: Vec<String> = if intervals.enabled {
        let mut v: Vec<String> = intervals.map.keys().cloned().collect();
        v.sort();
        v
    } else {
        let mut v: Vec<String> = fai.keys().cloned().collect();
        v.sort();
        v
    };
    contigs.dedup();

    let mut hist = [0u64; 101];
    let mut fa = File::open(&reference)?;
    for c in &contigs {
        if let Some(e) = fai.get(c) {
            let bases = load_contig(&mut fa, e)?;
            add_windows_by_gc(&bases, window, &mut hist);
        }
    }
    write_gc_index(&out_path, window, &ivhash, &hist)?;
    let total: u64 = hist.iter().sum();
    eprintln!("wrote {out_path}: window={window} intervals_sha={ivhash} total_windows={total}");
    Ok(())
}


/// Genome-wide per-base depth for the Lorenz-curve / Gini coefficient, reproducing
/// `bam-lorenz-coverage` (Hoogstrate et al., GigaScience 2021; GPL-3.0) which runs
/// `samtools depth` and builds a depth histogram over covered positions.
///
/// Matches `samtools depth` default counting: reads with UNMAP/SECONDARY/QCFAIL/DUP
/// are excluded (supplementary INCLUDED); depth is incremented over M/=/X reference
/// positions only (D/N advance without counting). No mapq/baseq filter, no depth cap.
///
/// Implemented as a delta-sweep per contig (coordinate-sorted BAM): O(reads) events
/// rather than a full per-base array, so no multi-GB buffer.
struct LorenzAcc {
    cur_ref: i64,
    deltas: std::collections::BTreeMap<i64, i64>, // pos -> net coverage delta (current contig)
    hist: Vec<u64>,                                // depth -> number of covered positions
    lengths: Vec<usize>,                           // contig lengths (by ref_id)
    seen: Vec<bool>,                               // contig received >=1 counted read
    investigated: u64,                             // sum of lengths of seen contigs (samtools depth -a domain)
}

impl LorenzAcc {
    fn new(lengths: Vec<usize>) -> Self {
        let n = lengths.len();
        LorenzAcc {
            cur_ref: -1,
            deltas: std::collections::BTreeMap::new(),
            hist: Vec::new(),
            lengths,
            seen: vec![false; n],
            investigated: 0,
        }
    }

    /// Sweep the current contig's delta events into the depth histogram.
    fn flush(&mut self) {
        if self.cur_ref < 0 {
            return;
        }
        let mut depth: i64 = 0;
        let mut last_pos: i64 = 0;
        for (&pos, &delta) in self.deltas.iter() {
            if depth > 0 && pos > last_pos {
                let len = (pos - last_pos) as u64;
                let d = depth as usize;
                if d >= self.hist.len() {
                    self.hist.resize(d + 1, 0);
                }
                self.hist[d] += len;
            }
            depth += delta;
            last_pos = pos;
        }
        self.deltas.clear();
        self.cur_ref = -1;
    }

    fn on_ref(&mut self, ref_id: i64) {
        if ref_id != self.cur_ref {
            self.flush();
            self.cur_ref = ref_id;
            if let Some(seen) = self.seen.get_mut(ref_id as usize) {
                if !*seen {
                    *seen = true;
                    self.investigated += *self.lengths.get(ref_id as usize).unwrap_or(&0) as u64;
                }
            }
        }
    }

    /// Add an aligned reference block [start0, end0) (0-based, end exclusive).
    fn add_block(&mut self, start0: i64, end0: i64) {
        let len = self.lengths.get(self.cur_ref as usize).copied().unwrap_or(0) as i64;
        let s = start0.max(0);
        let e = end0.min(len);
        if e > s {
            *self.deltas.entry(s).or_insert(0) += 1;
            *self.deltas.entry(e).or_insert(0) -= 1;
        }
    }

    /// Lorenz-curve ROC over covered positions (depth>=1), per bam-lorenz-coverage.
    /// Returns (roc, total_covered_positions, total_sequenced_bases). gini = 0.5 - roc.
    fn roc(&self) -> Option<(f64, u64, u64)> {
        let maxd = self.hist.len();
        let mut total_cov: u128 = 0;
        let mut total_bases: u128 = 0;
        for d in 1..maxd {
            let c = self.hist[d] as u128;
            total_cov += c;
            total_bases += c * d as u128;
        }
        if total_cov == 0 || total_bases == 0 {
            return None;
        }
        // Cumulative from highest depth to lowest, trapezoidal area (matches blc.py).
        let mut cumu_p: u128 = 0;
        let mut cumu_b: u128 = 0;
        let mut prev_p: u128 = 0;
        let mut prev_b: u128 = 0;
        let mut top: u128 = 0;
        for d in (1..maxd).rev() {
            if self.hist[d] == 0 {
                continue;
            }
            let c = self.hist[d] as u128;
            cumu_p += c;
            cumu_b += c * d as u128;
            let db = cumu_b - prev_b;
            top += db * (cumu_p - prev_p) + db * prev_p * 2;
            prev_p = cumu_p;
            prev_b = cumu_b;
        }
        let denom = total_bases * total_cov * 2;
        let roc = top as f64 / denom as f64;
        Some((roc, total_cov as u64, total_bases as u64))
    }
}


/// Optical-duplicate detector (Picard/Sentieon `READ_PAIR_OPTICAL_DUPLICATES`).
/// dupblaster marks duplicate *sets* but does not parse flow-cell coordinates, so we
/// recompute the optical subset here from read names. A duplicate set = read pairs sharing
/// the same 5' key (ref,pos,strand of read1 + mate ref/pos); within a set, read pairs whose
/// flow-cell (lane,tile,x,y) are within `pixel` of each other (transitively, same tile) are
/// optical. Count per set = sum over pixel-connected components of (component_size - 1),
/// matching Picard's OpticalDuplicateFinder. Reads are coordinate-sorted, so all members of
/// a set share read1's `pos` and are buffered together (bounded by pileup depth at a locus).
/// Read names without CASAVA `lane:tile:x:y` (non-Illumina / renamed FASTQ) contribute 0,
/// same as Picard.
struct OpticalAcc {
    pixel: i64,
    cur_ref: i64,
    cur_pos: i64,
    // read1-primary pairs at cur_pos: (mate_ref, mate_pos, reverse, lane, tile, x, y, has_xy)
    buf: Vec<(i64, i64, bool, u32, u32, i64, i64, bool)>,
    optical_pairs: u64,
}

impl OpticalAcc {
    fn new(pixel: i64) -> Self {
        OpticalAcc { pixel, cur_ref: -1, cur_pos: -1, buf: Vec::new(), optical_pairs: 0 }
    }

    /// Parse CASAVA 1.8 read name `<inst>:<run>:<fc>:<lane>:<tile>:<x>:<y>` -> (lane,tile,x,y).
    fn parse_xy(name: &[u8]) -> Option<(u32, u32, i64, i64)> {
        // read name is up to the first space; noodles already strips the comment.
        let s = std::str::from_utf8(name).ok()?;
        let f: Vec<&str> = s.split(':').collect();
        if f.len() < 7 { return None; }
        let n = f.len();
        let lane: u32 = f[n - 4].parse().ok()?;
        let tile: u32 = f[n - 3].parse().ok()?;
        let x: i64 = f[n - 2].parse().ok()?;
        // y can carry a trailing "/1" or "#index"; take leading integer
        let ypart = f[n - 1];
        let yend = ypart.find(|c: char| !c.is_ascii_digit()).unwrap_or(ypart.len());
        let y: i64 = ypart[..yend].parse().ok()?;
        Some((lane, tile, x, y))
    }

    fn flush(&mut self) {
        if self.buf.is_empty() { return; }
        // group by duplicate key tail (mate_ref, mate_pos, reverse)
        let mut items = std::mem::take(&mut self.buf);
        items.sort_unstable_by(|a, b| (a.0, a.1, a.2).cmp(&(b.0, b.1, b.2)));
        let mut i = 0;
        while i < items.len() {
            let mut j = i + 1;
            while j < items.len() && (items[j].0, items[j].1, items[j].2) == (items[i].0, items[i].1, items[i].2) {
                j += 1;
            }
            // set = items[i..j]; only members with flow-cell coords participate
            let set: Vec<&(i64, i64, bool, u32, u32, i64, i64, bool)> =
                items[i..j].iter().filter(|r| r.7).collect();
            let k = set.len();
            if k >= 2 {
                // union-find over same (lane,tile) within pixel distance
                let mut parent: Vec<usize> = (0..k).collect();
                fn find(p: &mut Vec<usize>, a: usize) -> usize {
                    let mut r = a; while p[r] != r { r = p[r]; }
                    let mut c = a; while p[c] != c { let n = p[c]; p[c] = r; c = n; } r
                }
                for a in 0..k {
                    for b in (a + 1)..k {
                        let (la, ta, xa, ya) = (set[a].3, set[a].4, set[a].5, set[a].6);
                        let (lb, tb, xb, yb) = (set[b].3, set[b].4, set[b].5, set[b].6);
                        if la == lb && ta == tb && (xa - xb).abs() <= self.pixel && (ya - yb).abs() <= self.pixel {
                            let ra = find(&mut parent, a); let rb = find(&mut parent, b);
                            if ra != rb { parent[ra] = rb; }
                        }
                    }
                }
                let mut comp_count = 0usize;
                for a in 0..k { if find(&mut parent, a) == a { comp_count += 1; } }
                self.optical_pairs += (k - comp_count) as u64; // Σ(component_size - 1)
            }
            i = j;
        }
    }

    fn add(&mut self, ref_id: i64, pos: i64, reverse: bool, mate_ref: i64, mate_pos: i64, name: &[u8]) {
        if ref_id != self.cur_ref || pos != self.cur_pos {
            self.flush();
            self.cur_ref = ref_id;
            self.cur_pos = pos;
        }
        match Self::parse_xy(name) {
            Some((lane, tile, x, y)) => self.buf.push((mate_ref, mate_pos, reverse, lane, tile, x, y, true)),
            None => self.buf.push((mate_ref, mate_pos, reverse, 0, 0, 0, 0, false)),
        }
    }

    fn merge(&mut self, o: &OpticalAcc) { self.optical_pairs += o.optical_pairs; }
}

/// Feed one record into an OpticalAcc: read1, primary, mapped, mate-mapped only.
fn optical_feed(record: &bam::Record, opt: &mut OpticalAcc) {
    let flags = record.flags();
    if flags.is_secondary() || flags.is_supplementary() || flags.is_unmapped()
        || !flags.is_segmented() || flags.is_mate_unmapped() || !flags.is_first_segment() {
        return;
    }
    let ref_id = record.reference_sequence_id().transpose().ok().flatten().map(|x| x as i64).unwrap_or(-1);
    let pos = record.alignment_start().and_then(|r| r.ok()).map(|p| usize::from(p) as i64).unwrap_or(-1);
    if ref_id < 0 || pos < 0 { return; }
    let mate_ref = record.mate_reference_sequence_id().transpose().ok().flatten().map(|x| x as i64).unwrap_or(-1);
    let mate_pos = record.mate_alignment_start().and_then(|r| r.ok()).map(|p| usize::from(p) as i64).unwrap_or(-1);
    let reverse = flags.is_reverse_complemented();
    let name = record.name();
    let name_bytes: &[u8] = name.as_ref().map(|n| n.as_ref()).unwrap_or(b"");
    opt.add(ref_id, pos, reverse, mate_ref, mate_pos, name_bytes);
}

/// Mergeable per-worker accumulator bundle for the region-sharded parallel path.
/// Every field is additive across disjoint contig sets, so N workers can each fill
/// their own `Acc` and the results merge into one identical to a single pass.
/// (Lorenz/preseq are intentionally excluded — they are low-pass-only and disabled
/// on the parallel WGS path.)
struct Acc {
    first: Cat,
    second: Cat,
    /// Single-end (non-segmented, FLAG 0x1 unset) reads, e.g. Ultima. Picard reports these
    /// under CATEGORY=UNPAIRED; they never enter the pair/chimera/insert accounting.
    unpaired: Cat,
    cov: CovAcc,
    insert_hist: Vec<u64>,
    insert_overflow: std::collections::BTreeMap<u64, u64>,
    insert_min: u64,
    insert_max: u64,
    chrm_aligned_bases: u64,
    qy_reads: u64,
    qy_total_bases: u64,
    qy_q20_bases: u64,
    qy_q30_bases: u64,
    /// Full base-quality histogram over the same bases QualityYield counts (index = phred,
    /// 0..=93 per the SAM spec). Q20/Q30 above are threshold counts and cannot be turned
    /// back into a distribution; HET_SNP_SENSITIVITY needs the whole shape, so keep it.
    /// Costs one increment in a loop that already reads every quality byte.
    qual_hist: [u64; QUAL_HIST_LEN],
    gcbias: Option<GcBias>,
    optical: OpticalAcc,
}

impl Acc {
    fn new(genome_territory: u64, coverage_cap: u32, gcbias: Option<GcBias>, optical_pixel: i64) -> Self {
        Acc {
            first: Cat::default(),
            second: Cat::default(),
            unpaired: Cat::default(),
            cov: CovAcc::new(genome_territory, coverage_cap, genome_territory > 0),
            insert_hist: vec![0u64; INSERT_CAP + 1],
            insert_overflow: std::collections::BTreeMap::new(),
            insert_min: u64::MAX,
            insert_max: 0,
            chrm_aligned_bases: 0,
            qy_reads: 0,
            qy_total_bases: 0,
            qy_q20_bases: 0,
            qy_q30_bases: 0,
            qual_hist: [0u64; QUAL_HIST_LEN],
            gcbias,
            optical: OpticalAcc::new(optical_pixel),
        }
    }

    fn merge(&mut self, mut o: Acc) {
        self.first.add(&o.first);
        self.second.add(&o.second);
        self.unpaired.add(&o.unpaired);
        self.cov.flush_all();
        o.cov.flush_all();
        self.cov.merge(&o.cov);
        for (i, c) in o.insert_hist.iter().enumerate() {
            self.insert_hist[i] += *c;
        }
        for (&value, &count) in o.insert_overflow.iter() {
            *self.insert_overflow.entry(value).or_insert(0) += count;
        }
        self.insert_min = self.insert_min.min(o.insert_min);
        self.insert_max = self.insert_max.max(o.insert_max);
        self.chrm_aligned_bases += o.chrm_aligned_bases;
        self.qy_reads += o.qy_reads;
        self.qy_total_bases += o.qy_total_bases;
        self.qy_q20_bases += o.qy_q20_bases;
        self.qy_q30_bases += o.qy_q30_bases;
        for i in 0..QUAL_HIST_LEN {
            self.qual_hist[i] += o.qual_hist[i];
        }
        if let (Some(dst), Some(src)) = (self.gcbias.as_mut(), o.gcbias.as_ref()) {
            dst.merge(src);
        }
        self.optical.flush();
        o.optical.flush();
        self.optical.merge(&o.optical);
    }
}

/// Per-record accumulation shared by the single-thread and parallel paths: QualityYield,
/// per-category AlignmentStat, WgsMetrics-style coverage, GC bias, chimera breakdown, and
/// insert size. Excludes Lorenz/preseq (low-pass only). Extracted verbatim from the
/// single-pass loop so the two paths produce identical numbers (validated on real BAMs).
fn process_record(record: &bam::Record, intervals: &Intervals, cov_intervals: &Intervals,
                  ref_names: &[String], acc: &mut Acc) {
    let Acc {
        first, second, unpaired, cov, insert_hist, insert_overflow, insert_min, insert_max, chrm_aligned_bases,
        qy_reads, qy_total_bases, qy_q20_bases, qy_q30_bases, qual_hist, gcbias, optical,
    } = acc;

    optical_feed(record, optical);

    let flags = record.flags();

    // primary reads only (skip secondary + supplementary)
    if flags.is_secondary() || flags.is_supplementary() {
        return;
    }

    // QualityYield: count over all primary reads (mapped or not), no interval filter,
    // matching Picard CollectQualityYieldMetrics. PF only (skip vendor-failed).
    if !flags.is_qc_fail() {
        *qy_reads += 1;
        let q = record.quality_scores();
        for &b in q.as_ref().iter() {
            *qy_total_bases += 1;
            if b >= 20 { *qy_q20_bases += 1; }
            if b >= 30 { *qy_q30_bases += 1; }
            // Clamp rather than index blindly: the SAM spec caps phred at 93, but a
            // malformed BAM must not panic a QC run.
            qual_hist[(b as usize).min(QUAL_HIST_LEN - 1)] += 1;
        }
    }

    // Paired reads need exactly one of READ1/READ2. Single-end reads (0x1 unset, e.g.
    // Ultima) go to the UNPAIRED category: same alignment/quality/coverage accounting,
    // no pair, chimera or insert-size contribution. Paired data never takes this branch,
    // so its output is unchanged.
    let is_paired = flags.is_segmented();
    let is_read1 = flags.is_first_segment();
    let is_read2 = flags.is_last_segment();
    if is_paired && !(is_read1 ^ is_read2) {
        return;
    }

    let mapped = !flags.is_unmapped();
    // 0x8 is only meaningful for paired reads; a single-end read has no mate to lose.
    let mate_unmapped = is_paired && flags.is_mate_unmapped();

    let (place_contig, place_pos): (Option<&str>, i64) = if mapped {
        let rid = record.reference_sequence_id().transpose().ok().flatten();
        let pos = record.alignment_start().and_then(|r| r.ok()).map(|p| usize::from(p) as i64).unwrap_or(0);
        (rid.and_then(|id| ref_names.get(id).map(|s| s.as_str())), pos)
    } else if is_paired && !mate_unmapped {
        let mrid = record.mate_reference_sequence_id().transpose().ok().flatten();
        let mpos = record.mate_alignment_start().and_then(|r| r.ok()).map(|p| usize::from(p) as i64).unwrap_or(0);
        (mrid.and_then(|id| ref_names.get(id).map(|s| s.as_str())), mpos)
    } else {
        (None, 0)
    };

    let in_interval = match place_contig {
        Some(c) => {
            let end1 = if mapped && intervals.overlap { alignment_end1(record, place_pos) } else { place_pos };
            intervals.admits(c, place_pos, end1)
        }
        None => false,
    };
    if !in_interval {
        return;
    }
    // Coverage is scoped separately (see --coverage-intervals): Sentieon's WgsMetricsAlgo
    // runs over a narrower region set than its AlignmentStat/GCBias pass, and matching that
    // split is the whole point of the second set. Safe to test only after the base-set gate
    // above because the coverage set is verified to be a subset of it at startup.
    let in_cov_interval = match place_contig {
        Some(c) => cov_intervals.contains(c, place_pos),
        None => false,
    };

    let cat = if !is_paired { unpaired } else if is_read1 { first } else { second };
    cat.total_reads += 1;
    let read_len = record.sequence().len() as u64;
    cat.read_length_sum += read_len;

    if mapped {
        cat.pf_reads_aligned += 1;
        if flags.is_reverse_complemented() {
            cat.neg_strand_aligned += 1;
        }
        let mapq = record.mapping_quality().map(u8::from).unwrap_or(0);
        let is_hq = mapq >= 20;
        let quals = record.quality_scores();
        let qual_bytes = quals.as_ref();
        let mut aln_bases: u64 = 0;
        let mut read_off: usize = 0;
        let mut hq_q20_aligned: u64 = 0;
        let mut ref_span: i64 = 0;
        {
            use noodles_sam::alignment::record::cigar::op::Kind;
            for op in record.cigar().iter().flatten() {
                let len = op.len();
                match op.kind() {
                    Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                        aln_bases += len as u64;
                        ref_span += len as i64;
                        if is_hq {
                            let end = (read_off + len).min(qual_bytes.len());
                            for &q in &qual_bytes[read_off.min(qual_bytes.len())..end] {
                                if q >= 20 { hq_q20_aligned += 1; }
                            }
                        }
                        read_off += len;
                    }
                    Kind::Insertion => {
                        cat.indel_events += 1;
                        cat.indel_bases += len as u64;
                        read_off += len;
                    }
                    Kind::SoftClip => { read_off += len; }
                    Kind::Deletion => {
                        cat.indel_events += 1;
                        cat.indel_bases += len as u64;
                        ref_span += len as i64;
                    }
                    Kind::Skip => { ref_span += len as i64; }
                    _ => {}
                }
            }
        }
        cat.aligned_bases += aln_bases;
        if let Some(c) = place_contig {
            if matches!(c, "chrM" | "MT" | "chrMT" | "M" | "chrM_rCRS") {
                *chrm_aligned_bases += aln_bases;
            }
        }

        if let Some(gb) = gcbias {
            if let Some(c) = place_contig {
                let start1 = place_pos;
                let end1 = place_pos + ref_span - 1;
                // a single-end read is its own cluster (Picard GcBias counts it as one)
                gb.add_read(c, start1, end1, flags.is_reverse_complemented(), is_read1 || !is_paired);
            }
        }

        let subs = if let Some(Ok(v)) = record.data().get(&Tag::from(*b"MD")) {
            match value_as_string(&v) {
                Some(md) => md_mismatches(&md),
                None => nm_minus_indels(record, aln_bases),
            }
        } else {
            nm_minus_indels(record, aln_bases)
        };
        cat.mismatch_sum += subs;

        if is_hq {
            cat.hq_reads += 1;
            cat.hq_aligned_bases += aln_bases;
            cat.hq_mismatch_sum += subs;
            cat.hq_q20_bases += hq_q20_aligned;
        }
        // Single-end chimeras (Picard fragment rule): a MAPQ>=20 mapped read is chimeric
        // when it carries a supplementary alignment (SA tag).
        if !is_paired && is_hq {
            cat.chim_fragment_den += 1;
            if matches!(record.data().get(&Tag::from(*b"SA")), Some(Ok(_))) {
                cat.chim_split_sa += 1;
                cat.chimeric_reads += 1;
            }
        }

        if cov.enabled && in_cov_interval {
            let ref_id0 = record.reference_sequence_id().transpose().ok().flatten().map(|x| x as i64).unwrap_or(-1);
            let start0 = record.alignment_start().and_then(|r| r.ok()).map(|p| usize::from(p) as i64 - 1).unwrap_or(0);
            let is_dup = (u16::from(flags) & 0x400) != 0;
            if !is_hq {
                cov.exc_mapq += aln_bases;
            } else if is_dup {
                cov.exc_dupe += aln_bases;
            } else if mate_unmapped {
                cov.exc_unpaired += aln_bases;
            } else {
                cov.on_new_read(ref_id0, start0);
                let mate_start0 = record.mate_alignment_start().and_then(|r| r.ok()).map(|p| usize::from(p) as i64 - 1).unwrap_or(-1);
                let mate_ref = record.mate_reference_sequence_id().transpose().ok().flatten().map(|x| x as i64).unwrap_or(-2);
                let is_later = start0 > mate_start0 || (start0 == mate_start0 && is_read2);
                let (ov_s, ov_e) = if mate_ref == ref_id0 && mate_start0 >= 0 && is_later {
                    let mc_len = match record.data().get(&Tag::from(*b"MC")) {
                        Some(Ok(v)) => value_as_string(&v).map(|s| ref_len_from_cigar_str(&s)).unwrap_or(0),
                        _ => 0,
                    };
                    (mate_start0, mate_start0 + mc_len)
                } else {
                    (0, -1)
                };
                let mut rpos = start0;
                let mut roff: usize = 0;
                use noodles_sam::alignment::record::cigar::op::Kind;
                for op in record.cigar().iter().flatten() {
                    let len = op.len();
                    match op.kind() {
                        Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                            for k in 0..len {
                                let q = qual_bytes.get(roff + k).copied().unwrap_or(0);
                                let p = rpos + k as i64;
                                if q < 20 {
                                    cov.exc_baseq += 1;
                                } else if p >= ov_s && p < ov_e {
                                    cov.exc_overlap += 1;
                                } else {
                                    cov.add_position(p);
                                }
                            }
                            rpos += len as i64;
                            roff += len;
                        }
                        Kind::Insertion | Kind::SoftClip => { roff += len; }
                        Kind::Deletion | Kind::Skip => { rpos += len as i64; }
                        _ => {}
                    }
                }
            }
        }

        if is_paired && !mate_unmapped {
            cat.reads_aligned_in_pairs += 1;
            let rid = record.reference_sequence_id().transpose().ok().flatten();
            let mrid = record.mate_reference_sequence_id().transpose().ok().flatten();
            let same_contig = matches!((rid, mrid), (Some(a), Some(b)) if a == b);
            let tlen = i32::from(record.template_length()) as i64;
            let read_rev = flags.is_reverse_complemented();
            let mate_rev = flags.is_mate_reverse_complemented();
            let fr_proper = (!read_rev && mate_rev && tlen >= 0) || (read_rev && !mate_rev && tlen <= 0);
            let within_insert = tlen.abs() <= 100_000;
            let has_sa = matches!(record.data().get(&Tag::from(*b"SA")), Some(Ok(_)));

            let diff_contig = !same_contig;
            let bad_orient = same_contig && !fr_proper;
            let large_insert = same_contig && fr_proper && !within_insert;
            let chimeric = diff_contig || bad_orient || large_insert || has_sa;

            if diff_contig { cat.chim_diff_contig += 1; }
            if bad_orient { cat.chim_bad_orient += 1; classify_bad_orient(cat, read_rev, mate_rev, tlen); }
            if large_insert { cat.chim_large_insert += 1; }
            if has_sa { cat.chim_split_sa += 1; }
            if chimeric { cat.chimeric_reads += 1; }
            cat.chim_cube.0[chim_cube_bin(same_contig, fr_proper, read_rev, mate_rev, tlen, has_sa)] += 1;

            let properly_paired = (u16::from(flags) & 0x2) != 0;
            if properly_paired && fr_proper && tlen > 0 {
                let insert = tlen as u64;
                *insert_min = (*insert_min).min(insert);
                *insert_max = (*insert_max).max(insert);
                if insert <= INSERT_CAP as u64 {
                    insert_hist[insert as usize] += 1;
                } else {
                    *insert_overflow.entry(insert).or_insert(0) += 1;
                }
            }
        }
    }
}

/// Region-sharded parallel accumulation (WGS path). Each worker pulls whole contigs
/// (largest-first for load balance) from a shared queue, queries them from its own
/// indexed reader, and accumulates into a private `Acc`; one work unit processes the
/// unmapped reads at the end of the file (needed for QualityYield totals). All `Acc`s
/// merge into one identical to a single sequential pass. Coordinate order is preserved
/// within each contig, so the coverage pileup is correct per shard.
///
/// Memory is bounded: per-worker footprint = coverage histogram + insert histogram
/// (~8 MB) + at most one contig's GC-bias base buffer (only when --reference is given).
/// Total is therefore ~workers x (few MB + one contig), independent of read depth.
fn run_parallel_regions(
    bam_path: &str,
    intervals: &Intervals,
    cov_intervals: &Intervals,
    ref_names: &[String],
    ref_lengths: &[usize],
    genome_territory: u64,
    coverage_cap: u32,
    gc_template: Option<&GcBias>,
    threads: usize,
    optical_pixel: i64,
) -> std::io::Result<Acc> {
    use std::sync::atomic::{AtomicUsize, Ordering};
    use noodles_core::Region;

    let n_contigs = ref_names.len();
    // Work order: contig indices sorted by length descending, then a sentinel (== n_contigs)
    // for the unmapped-reads unit.
    let mut order: Vec<usize> = (0..n_contigs).collect();
    order.sort_by(|&a, &b| ref_lengths[b].cmp(&ref_lengths[a]));
    order.push(n_contigs); // unmapped sentinel

    let next = AtomicUsize::new(0);
    let n_workers = threads.max(1).min(order.len());

    let build_gc = |contigs: &[String]| -> Option<GcBias> {
        // windows_cache is Some, so this only reloads the .fai (already validated when
        // the template was built) — no full-genome scan.
        gc_template.and_then(|t| GcBias::new(&t.fa_path, t.window, contigs, Some(t.windows_by_gc)).ok())
    };

    let merged = std::thread::scope(|scope| -> std::io::Result<Acc> {
        let mut handles = Vec::new();
        for _ in 0..n_workers {
            let next = &next;
            let order = &order;
            let build_gc = &build_gc;
            handles.push(scope.spawn(move || -> std::io::Result<Acc> {
                let mut acc = Acc::new(genome_territory, coverage_cap, build_gc(ref_names), optical_pixel);
                loop {
                    let slot = next.fetch_add(1, Ordering::Relaxed);
                    if slot >= order.len() { break; }
                    let unit = order[slot];
                    let mut reader = bam::io::indexed_reader::Builder::default()
                        .build_from_path(bam_path)?;
                    let header = reader.read_header()?;
                    if unit == n_contigs {
                        // unmapped reads (no position) — only contribute QualityYield.
                        for result in reader.query_unmapped()? {
                            let record = result?;
                            process_record(&record, intervals, cov_intervals, ref_names, &mut acc);
                        }
                    } else {
                        // Construct the region directly (do NOT string-parse: contig
                        // names can contain ':' / '*', e.g. HLA/alt contigs, which the
                        // "name:start-end" grammar would misinterpret). Whole contig.
                        let region = Region::new(ref_names[unit].as_bytes(), ..);
                        let query = reader.query(&header, &region)?;
                        for result in query {
                            let record = result?;
                            process_record(&record, intervals, cov_intervals, ref_names, &mut acc);
                        }
                    }
                    acc.cov.flush_all();
                }
                Ok(acc)
            }));
        }
        let mut merged = Acc::new(genome_territory, coverage_cap, build_gc(ref_names), optical_pixel);
        for h in handles {
            match h.join() {
                Ok(Ok(acc)) => merged.merge(acc),
                Ok(Err(e)) => return Err(e),
                Err(_) => return Err(std::io::Error::new(std::io::ErrorKind::Other, "worker panicked")),
            }
        }
        Ok(merged)
    })?;

    Ok(merged)
}

fn run() -> std::io::Result<()> {
    // ---- args ----
    let args: Vec<String> = std::env::args().collect();
    let mut bam_path = String::new();
    let mut bed_path: Option<String> = None;
    // Region set for the WgsMetrics-style coverage block only. Sentieon runs AlignmentStat /
    // GCBias / InsertSize over base_metrics_intervals (chr1-22 + X + Y + M) but WgsMetricsAlgo
    // over the narrower wgs_or_target_intervals (chr1-22), and genome_territory is the non-N
    // size of that narrower set. Defaults to --intervals, reproducing the single-set behaviour.
    let mut cov_bed_path: Option<String> = None;
    let mut sample = String::from("sample");
    let mut out_path = String::from("-");
    let mut genome_territory: u64 = 0;
    let mut coverage_cap: u32 = 250;
    let mut reference: Option<String> = None;
    let mut gc_windows: Option<String> = None;
    let mut threads: usize = 1;
    // preseq (library complexity) and Lorenz/Gini (coverage evenness) are low-pass
    // single-cell QC concepts (used by basej-dnaqc on 2M-read subsamples). They are
    // meaningless on high-pass WGS and dominate runtime at full depth, so allow the
    // WGS path to disable them. Default ON to keep basej-dnaqc-rs behavior unchanged.
    let mut run_preseq = true;
    let mut run_lorenz = true;
    let mut interval_overlap = false;
    // Optical-duplicate pixel distance (Picard OpticalDuplicateFinder default 100 for
    // unpatterned flowcells; use 2500 for patterned NovaSeq/HiSeqX). Read pairs in the
    // same duplicate set whose flow-cell (x,y) are within this many pixels are optical.
    let mut optical_pixel: i64 = 100;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--bam" => { bam_path = args[i + 1].clone(); i += 2; }
            "--intervals" => { bed_path = Some(args[i + 1].clone()); i += 2; }
            "--coverage-intervals" => { cov_bed_path = Some(args[i + 1].clone()); i += 2; }
            "--sample" => { sample = args[i + 1].clone(); i += 2; }
            "--out" => { out_path = args[i + 1].clone(); i += 2; }
            "--genome-territory" => { genome_territory = args[i + 1].parse().unwrap_or(0); i += 2; }
            "--coverage-cap" => { coverage_cap = args[i + 1].parse().unwrap_or(250); i += 2; }
            "--reference" => { reference = Some(args[i + 1].clone()); i += 2; }
            "--gc-windows" => { gc_windows = Some(args[i + 1].clone()); i += 2; }
            "--threads" | "-t" => { threads = args[i + 1].parse().unwrap_or(1).max(1); i += 2; }
            "--no-preseq" => { run_preseq = false; i += 1; }
            "--no-lorenz" => { run_lorenz = false; i += 1; }
            "--interval-overlap" => { interval_overlap = true; i += 1; }
            "--no-lowpass" => { run_preseq = false; run_lorenz = false; i += 1; }
            "--optical-pixel-distance" => { optical_pixel = args[i + 1].parse().unwrap_or(100); i += 2; }
            "-h" | "--help" => {
                eprintln!("usage: bskryb-qc --bam <in.bam> [--intervals <in.bed>] [--coverage-intervals <cov.bed>] [--sample <name>] [--out <out.json|->] [--genome-territory <N>] [--coverage-cap <N>] [--reference <genome.fa>] [--gc-windows <index.txt>] [--threads <N>] [--no-preseq] [--no-lorenz] [--no-lowpass] [--optical-pixel-distance <N>] [--interval-overlap]");
                return Ok(());
            }
            other => { eprintln!("unknown arg: {other}"); std::process::exit(2); }
        }
    }
    if bam_path.is_empty() {
        eprintln!("--bam is required");
        std::process::exit(2);
    }

    let mut intervals = match &bed_path {
        Some(p) => Intervals::from_bed(p)?,
        None => Intervals::none(),
    };
    // Only the base set's placement rule changes; coverage stays per-base.
    intervals.overlap = interval_overlap;
    // Omitting --coverage-intervals reuses the base set, i.e. the pre-split behaviour.
    let cov_intervals = match &cov_bed_path {
        Some(p) => Intervals::from_bed(p)?,
        None => intervals.clone(),
    };
    if let Some(why) = coverage_subset_violation(&cov_intervals, &intervals) {
        eprintln!("bskryb-qc error: {why}");
        eprintln!(
            "  --coverage-intervals must be a subset of --intervals; coverage is evaluated \
             only for reads that already passed the --intervals filter."
        );
        std::process::exit(2);
    }

    // ---- open BAM ----
    // Multithreaded BGZF decompression when --threads > 1 (decompression is the
    // dominant cost of this counting pass); single-threaded reader otherwise.
    let mut reader = {
        use std::num::NonZeroUsize;
        let file = File::open(&bam_path)?;
        let worker_count = NonZeroUsize::new(threads).unwrap_or(NonZeroUsize::MIN);
        let decoder = noodles_bgzf::MultithreadedReader::with_worker_count(worker_count, file);
        bam::io::Reader::from(decoder)
    };
    let header = reader.read_header()?;

    // ref_id -> contig name
    let ref_names: Vec<String> = header
        .reference_sequences()
        .keys()
        .map(|k| k.to_string())
        .collect();

    // ref_id -> contig length (for the Lorenz/Gini genome-wide depth accumulator)
    let ref_lengths: Vec<usize> = header
        .reference_sequences()
        .values()
        .map(|rs| usize::from(rs.length()))
        .collect();
    let mut lorenz = LorenzAcc::new(ref_lengths.clone());

    let mut first = Cat::default();
    let mut second = Cat::default();
    let mut unpaired = Cat::default();
    // Insert sizes use a dense histogram through INSERT_CAP plus a sparse exact
    // overflow map. This keeps memory independent of pair count while preserving
    // canonical median/MAD/min/max even for unusually large valid TLEN values.
    let mut insert_hist: Vec<u64> = vec![0u64; INSERT_CAP + 1]; // exact TLEN <= cap
    let mut insert_overflow: std::collections::BTreeMap<u64, u64> =
        std::collections::BTreeMap::new();
    let mut insert_min: u64 = u64::MAX;
    let mut insert_max: u64 = 0;
    let mut cov = CovAcc::new(genome_territory, coverage_cap, genome_territory > 0);
    let mut chrm_aligned_bases: u64 = 0; // aligned bases on chrM/MT (for pct_chrm)
    let mut preseq_acc = PreseqAcc::new(10); // gc_extrap coverage histogram (10bp bins)

    // GC bias (Picard CollectGcBiasMetrics-equivalent), enabled when --reference given.
    // Uses a cached windowsByGc index (--gc-windows) when provided (single-pass, no
    // full-genome scan); otherwise computes windowsByGc over the interval contigs.
    let mut gcbias = match &reference {
        Some(fa) => {
            let contigs: Vec<String> = if intervals.enabled {
                let mut v: Vec<String> = intervals.map.keys().cloned().collect();
                v.sort();
                v
            } else {
                ref_names.clone()
            };
            // Load cache and validate it matches the window size + interval set in use.
            let windows_cache = match &gc_windows {
                Some(idx) => match read_gc_index(idx) {
                    Ok((hist, win, ivhash)) => {
                        let cur_ivhash = if intervals.enabled { intervals_hash(&intervals) } else { "none".to_string() };
                        if win != 100 {
                            eprintln!("gc-windows cache window={win} != 100; recomputing from reference");
                            None
                        } else if ivhash != cur_ivhash {
                            eprintln!("gc-windows cache intervals_sha mismatch ({ivhash} != {cur_ivhash}); recomputing from reference");
                            None
                        } else {
                            Some(hist)
                        }
                    }
                    Err(e) => { eprintln!("could not read gc-windows cache ({e}); recomputing"); None }
                },
                None => None,
            };
            match GcBias::new(fa, 100, &contigs, windows_cache) {
                Ok(g) => Some(g),
                Err(e) => { eprintln!("gcbias disabled: {e}"); None }
            }
        }
        None => None,
    };
    // QualityYield (Picard CollectQualityYieldMetrics-equivalent): over all primary reads
    let mut qy_reads: u64 = 0;
    let mut qy_total_bases: u64 = 0;
    let mut qy_q20_bases: u64 = 0;
    let mut qy_q30_bases: u64 = 0;
    let mut qual_hist: [u64; QUAL_HIST_LEN] = [0u64; QUAL_HIST_LEN];

    // Optical-duplicate accumulator (Picard/Sentieon READ_PAIR_OPTICAL_DUPLICATES).
    let mut optical = OpticalAcc::new(optical_pixel);

    // Parallel region-sharded path for high-pass WGS (only when the low-pass metrics
    // are disabled). The validated single-thread loop below handles everything else,
    // including all basej-dnaqc-rs (low-pass) runs, and is left byte-for-byte unchanged.
    let parallel = threads > 1 && !run_preseq && !run_lorenz;
    if parallel {
        let mut merged = run_parallel_regions(
            &bam_path, &intervals, &cov_intervals, &ref_names, &ref_lengths,
            genome_territory, coverage_cap, gcbias.as_ref(), threads, optical_pixel,
        )?;
        first = merged.first;
        second = merged.second;
        unpaired = merged.unpaired;
        cov = merged.cov;
        insert_hist = merged.insert_hist;
        insert_overflow = merged.insert_overflow;
        insert_min = merged.insert_min;
        insert_max = merged.insert_max;
        chrm_aligned_bases = merged.chrm_aligned_bases;
        qy_reads = merged.qy_reads;
        qy_total_bases = merged.qy_total_bases;
        qy_q20_bases = merged.qy_q20_bases;
        qy_q30_bases = merged.qy_q30_bases;
        qual_hist = merged.qual_hist;
        gcbias = merged.gcbias;
        merged.optical.flush();
        optical.merge(&merged.optical);
    } else {
    for result in reader.records() {
        let record = result?;
        let flags = record.flags();

        // ---- optical-duplicate detection (read1/primary/mapped/mate-mapped only) ----
        optical_feed(&record, &mut optical);

        // ---- Lorenz/Gini genome-wide depth (samtools depth default semantics) ----
        // Include supplementary; exclude UNMAP/SECONDARY/QCFAIL/DUP. Must run BEFORE
        // the primary-only skip below (supplementary alignments contribute to depth).
        // Skipped on the high-pass WGS path (--no-lorenz): evenness is a low-pass metric.
        if run_lorenz && !(flags.is_unmapped() || flags.is_secondary() || flags.is_qc_fail() || flags.is_duplicate()) {
            if let Some(ref_id) = record.reference_sequence_id().transpose().ok().flatten() {
                lorenz.on_ref(ref_id as i64);
                if let Some(Ok(start)) = record.alignment_start() {
                    let mut pos = usize::from(start) as i64 - 1; // 0-based
                    use noodles_sam::alignment::record::cigar::op::Kind;
                    for op in record.cigar().iter().flatten() {
                        let l = op.len() as i64;
                        match op.kind() {
                            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                                lorenz.add_block(pos, pos + l);
                                pos += l;
                            }
                            Kind::Deletion | Kind::Skip => { pos += l; }
                            _ => {}
                        }
                    }
                }
            }
        }

        // primary reads only (skip secondary + supplementary)
        if flags.is_secondary() || flags.is_supplementary() {
            continue;
        }

        // QualityYield: count over all primary reads (mapped or not), no interval filter,
        // matching Picard CollectQualityYieldMetrics. PF only (skip vendor-failed).
        if !flags.is_qc_fail() {
            qy_reads += 1;
            let q = record.quality_scores();
            for &b in q.as_ref().iter() {
                qy_total_bases += 1;
                if b >= 20 { qy_q20_bases += 1; }
                if b >= 30 { qy_q30_bases += 1; }
                qual_hist[(b as usize).min(QUAL_HIST_LEN - 1)] += 1;
            }
        }

        // preseq gc_extrap coverage binning: all primary mapped reads, no interval filter
        // (matches bam2mr -> gc_extrap, which is not interval-restricted).
        // preseq library-complexity binning (low-pass metric; skipped with --no-preseq).
        if run_preseq && !flags.is_unmapped() {
            let ref_id0 = record
                .reference_sequence_id()
                .transpose()
                .ok()
                .flatten()
                .map(|x| x as i64)
                .unwrap_or(-1);
            let start1 = record
                .alignment_start()
                .and_then(|r| r.ok())
                .map(|p| usize::from(p) as i64)
                .unwrap_or(0);
            let mut rspan: i64 = 0;
            {
                use noodles_sam::alignment::record::cigar::op::Kind;
                for op in record.cigar().iter().flatten() {
                    match op.kind() {
                        Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch
                        | Kind::Deletion | Kind::Skip => rspan += op.len() as i64,
                        _ => {}
                    }
                }
            }
            // mate-overlap region to exclude (fragment merge): only for the later mate
            // of a proper pair on the same contig, using the MC tag for the mate's end.
            let mut ov_s: i64 = 0;
            let mut ov_e: i64 = -1;
            let properly_paired = (u16::from(flags) & 0x2) != 0;
            if properly_paired && !flags.is_mate_unmapped() {
                let mate_ref = record.mate_reference_sequence_id().transpose().ok().flatten().map(|x| x as i64).unwrap_or(-2);
                if mate_ref == ref_id0 {
                    let mate_start0 = record
                        .mate_alignment_start()
                        .and_then(|r| r.ok())
                        .map(|p| usize::from(p) as i64 - 1)
                        .unwrap_or(-1);
                    let start0 = start1 - 1;
                    let is_read2 = flags.is_last_segment();
                    let is_later = start0 > mate_start0 || (start0 == mate_start0 && is_read2);
                    if mate_start0 >= 0 && is_later {
                        let mc_len = match record.data().get(&Tag::from(*b"MC")) {
                            Some(Ok(v)) => value_as_string(&v).map(|s| ref_len_from_cigar_str(&s)).unwrap_or(0),
                            _ => 0,
                        };
                        if mc_len > 0 {
                            ov_s = mate_start0;
                            ov_e = mate_start0 + mc_len;
                        }
                    }
                }
            }
            preseq_acc.add_read(ref_id0, start1, rspan, ov_s, ov_e);
        }

        // Paired reads need a defined mate segment role; single-end reads (0x1 unset,
        // e.g. Ultima) are accumulated under UNPAIRED (see process_record).
        let is_paired = flags.is_segmented();
        let is_read1 = flags.is_first_segment();
        let is_read2 = flags.is_last_segment();
        if is_paired && !(is_read1 ^ is_read2) {
            continue;
        }

        let mapped = !flags.is_unmapped();
        let mate_unmapped = is_paired && flags.is_mate_unmapped();

        // Determine placement contig+pos for the interval test.
        // Mapped read -> own position. Placed-unmapped read -> mate position.
        let (place_contig, place_pos): (Option<&str>, i64) = if mapped {
            let rid = record.reference_sequence_id().transpose().ok().flatten();
            let pos = record
                .alignment_start()
                .and_then(|r| r.ok())
                .map(|p| usize::from(p) as i64)
                .unwrap_or(0);
            (rid.and_then(|id| ref_names.get(id).map(|s| s.as_str())), pos)
        } else if is_paired && !mate_unmapped {
            let mrid = record.mate_reference_sequence_id().transpose().ok().flatten();
            let mpos = record
                .mate_alignment_start()
                .and_then(|r| r.ok())
                .map(|p| usize::from(p) as i64)
                .unwrap_or(0);
            (mrid.and_then(|id| ref_names.get(id).map(|s| s.as_str())), mpos)
        } else {
            (None, 0)
        };

        let in_interval = match place_contig {
            Some(c) => {
                let end1 = if mapped && intervals.overlap { alignment_end1(&record, place_pos) } else { place_pos };
                intervals.admits(c, place_pos, end1)
            }
            None => false,
        };
        if !in_interval {
            continue;
        }
        // See process_record: coverage is scoped by its own (subset) interval list.
        let in_cov_interval = match place_contig {
            Some(c) => cov_intervals.contains(c, place_pos),
            None => false,
        };

        let cat = if !is_paired { &mut unpaired } else if is_read1 { &mut first } else { &mut second };
        cat.total_reads += 1;
        let read_len = record.sequence().len() as u64;
        cat.read_length_sum += read_len;

        if mapped {
            cat.pf_reads_aligned += 1;
            if flags.is_reverse_complemented() {
                cat.neg_strand_aligned += 1;
            }
            // Single CIGAR walk: aligned bases (M/=/X), indel events/bases, and
            // (for HQ reads) Q20 bases restricted to aligned read positions.
            let properly_paired = (u16::from(flags) & 0x2) != 0;
            let mapq = record.mapping_quality().map(u8::from).unwrap_or(0);
            let is_hq = mapq >= 20;
            let quals = record.quality_scores();
            let qual_bytes = quals.as_ref();
            let mut aln_bases: u64 = 0;
            let mut read_off: usize = 0;
            let mut hq_q20_aligned: u64 = 0;
            let mut ref_span: i64 = 0;
            {
                use noodles_sam::alignment::record::cigar::op::Kind;
                for op in record.cigar().iter().flatten() {
                    let len = op.len();
                    match op.kind() {
                        Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                            aln_bases += len as u64;
                            ref_span += len as i64;
                            if is_hq {
                                let end = (read_off + len).min(qual_bytes.len());
                                for &q in &qual_bytes[read_off.min(qual_bytes.len())..end] {
                                    if q >= 20 { hq_q20_aligned += 1; }
                                }
                            }
                            read_off += len;
                        }
                        Kind::Insertion => {
                            cat.indel_events += 1;
                            cat.indel_bases += len as u64;
                            read_off += len;
                        }
                        Kind::SoftClip => {
                            read_off += len;
                        }
                        Kind::Deletion => {
                            cat.indel_events += 1;
                            cat.indel_bases += len as u64;
                            ref_span += len as i64;
                        }
                        Kind::Skip => {
                            ref_span += len as i64;
                        }
                        _ => {}
                    }
                }
            }
            cat.aligned_bases += aln_bases;
            if let Some(c) = place_contig {
                if matches!(c, "chrM" | "MT" | "chrMT" | "M" | "chrM_rCRS") {
                    chrm_aligned_bases += aln_bases;
                }
            }

            // GC bias: assign read to the window at its 5' start (Picard convention).
            if let Some(ref mut gb) = gcbias {
                if let Some(c) = place_contig {
                    let start1 = place_pos; // 1-based alignment start
                    let end1 = place_pos + ref_span - 1; // 1-based inclusive end
                    gb.add_read(c, start1, end1, flags.is_reverse_complemented(), is_read1 || !is_paired);
                }
            }

            // substitutions: MD-tag substitution count (== NM - indel_bases when consistent)
            let subs = if let Some(Ok(v)) = record.data().get(&Tag::from(*b"MD")) {
                match value_as_string(&v) {
                    Some(md) => md_mismatches(&md),
                    None => nm_minus_indels(&record, aln_bases),
                }
            } else {
                nm_minus_indels(&record, aln_bases)
            };
            cat.mismatch_sum += subs;

            if is_hq {
                cat.hq_reads += 1;
                cat.hq_aligned_bases += aln_bases;
                cat.hq_mismatch_sum += subs;
                cat.hq_q20_bases += hq_q20_aligned;
            }
            // Single-end chimeras (Picard fragment rule): a MAPQ>=20 mapped read is chimeric
            // when it carries a supplementary alignment (SA tag).
            if !is_paired && is_hq {
                cat.chim_fragment_den += 1;
                if matches!(record.data().get(&Tag::from(*b"SA")), Some(Ok(_))) {
                    cat.chim_split_sa += 1;
                    cat.chimeric_reads += 1;
                }
            }

            // ---- WgsMetrics-style coverage (single-pass pileup) ----
            if cov.enabled && in_cov_interval {
                let ref_id0 = record
                    .reference_sequence_id()
                    .transpose()
                    .ok()
                    .flatten()
                    .map(|x| x as i64)
                    .unwrap_or(-1);
                let start0 = record
                    .alignment_start()
                    .and_then(|r| r.ok())
                    .map(|p| usize::from(p) as i64 - 1)
                    .unwrap_or(0);
                let is_dup = (u16::from(flags) & 0x400) != 0;
                if !is_hq {
                    cov.exc_mapq += aln_bases;
                } else if is_dup {
                    cov.exc_dupe += aln_bases;
                } else if mate_unmapped {
                    cov.exc_unpaired += aln_bases;
                } else {
                    cov.on_new_read(ref_id0, start0);
                    // mate overlap region (clip the later-starting mate) via MC tag
                    let mate_start0 = record
                        .mate_alignment_start()
                        .and_then(|r| r.ok())
                        .map(|p| usize::from(p) as i64 - 1)
                        .unwrap_or(-1);
                    let mate_ref = record
                        .mate_reference_sequence_id()
                        .transpose()
                        .ok()
                        .flatten()
                        .map(|x| x as i64)
                        .unwrap_or(-2);
                    let is_later = start0 > mate_start0 || (start0 == mate_start0 && is_read2);
                    let (ov_s, ov_e) = if mate_ref == ref_id0 && mate_start0 >= 0 && is_later {
                        let mc_len = match record.data().get(&Tag::from(*b"MC")) {
                            Some(Ok(v)) => value_as_string(&v).map(|s| ref_len_from_cigar_str(&s)).unwrap_or(0),
                            _ => 0,
                        };
                        (mate_start0, mate_start0 + mc_len)
                    } else {
                        (0, -1)
                    };
                    // walk CIGAR tracking reference position and per-base quality
                    let mut rpos = start0;
                    let mut roff: usize = 0;
                    use noodles_sam::alignment::record::cigar::op::Kind;
                    for op in record.cigar().iter().flatten() {
                        let len = op.len();
                        match op.kind() {
                            Kind::Match | Kind::SequenceMatch | Kind::SequenceMismatch => {
                                for k in 0..len {
                                    let q = qual_bytes.get(roff + k).copied().unwrap_or(0);
                                    let p = rpos + k as i64;
                                    if q < 20 {
                                        cov.exc_baseq += 1;
                                    } else if p >= ov_s && p < ov_e {
                                        cov.exc_overlap += 1;
                                    } else {
                                        cov.add_position(p);
                                    }
                                }
                                rpos += len as i64;
                                roff += len;
                            }
                            Kind::Insertion | Kind::SoftClip => {
                                roff += len;
                            }
                            Kind::Deletion | Kind::Skip => {
                                rpos += len as i64;
                            }
                            _ => {}
                        }
                    }
                }
            }

            // Picard chimera + component breakdown (paired reads only)
            if is_paired && !mate_unmapped {
                cat.reads_aligned_in_pairs += 1;
                let rid = record.reference_sequence_id().transpose().ok().flatten();
                let mrid = record.mate_reference_sequence_id().transpose().ok().flatten();
                let same_contig = matches!((rid, mrid), (Some(a), Some(b)) if a == b);
                let tlen = i32::from(record.template_length()) as i64;
                let read_rev = flags.is_reverse_complemented();
                let mate_rev = flags.is_mate_reverse_complemented();
                let fr_proper = (!read_rev && mate_rev && tlen >= 0)
                    || (read_rev && !mate_rev && tlen <= 0);
                let within_insert = tlen.abs() <= 100_000;
                let has_sa = matches!(record.data().get(&Tag::from(*b"SA")), Some(Ok(_)));

                let diff_contig = !same_contig;
                let bad_orient = same_contig && !fr_proper;
                let large_insert = same_contig && fr_proper && !within_insert;
                let chimeric = diff_contig || bad_orient || large_insert || has_sa;

                if diff_contig { cat.chim_diff_contig += 1; }
                if bad_orient { cat.chim_bad_orient += 1; classify_bad_orient(cat, read_rev, mate_rev, tlen); }
                if large_insert { cat.chim_large_insert += 1; }
                if has_sa { cat.chim_split_sa += 1; }
                if chimeric { cat.chimeric_reads += 1; }
                cat.chim_cube.0[chim_cube_bin(same_contig, fr_proper, read_rev, mate_rev, tlen, has_sa)] += 1;

                // insert size: count each proper FR pair once via its forward mate (TLEN>0)
                if properly_paired && fr_proper && tlen > 0 {
                    let insert = tlen as u64;
                    insert_min = insert_min.min(insert);
                    insert_max = insert_max.max(insert);
                    if insert <= INSERT_CAP as u64 {
                        insert_hist[insert as usize] += 1;
                    } else {
                        *insert_overflow.entry(insert).or_insert(0) += 1;
                    }
                }
            }
        }
    }
    } // end else: single-thread record loop

    let mut pair = first.clone();
    pair.add(&second);
    cov.flush_all();
    optical.flush(); // finalize any buffered duplicate set (single-thread path)
    let optical_dups = optical.optical_pairs;
    // Lorenz/Gini (low-pass evenness) — only finalize when enabled.
    let lorenz_roc = if run_lorenz {
        lorenz.flush();
        lorenz.roc() // (roc, total_covered, total_bases); gini = 0.5 - roc
    } else {
        None
    };
    // preseq (low-pass library complexity) — only finalize when enabled.
    let (preseq_total_bins, preseq_distinct_bins, preseq_max_cov, preseq_count) = if run_preseq {
        preseq_acc.flush_all();
        let tb: f64 = preseq_acc.hist.iter().enumerate().map(|(i, &c)| i as f64 * c).sum();
        let db: f64 = preseq_acc.hist.iter().sum();
        let mc = preseq_acc.hist.len().saturating_sub(1);
        let pc = preseq::gc_extrap_estimate(&preseq_acc.hist, 10, 1.0e8, 1.0e12, 100, 408);
        (tb, db, mc, pc)
    } else {
        (0.0, 0.0, 0usize, None)
    };

    // ---- insert size stats (FR) ----
    let (ins_median, ins_mad, ins_mean, ins_sd, ins_pairs) =
        insert_stats_hist(&insert_hist, &insert_overflow);
    let ins_min = if ins_pairs > 0 { insert_min } else { 0 };
    let ins_max = if ins_pairs > 0 { insert_max } else { 0 };

    // ---- emit JSON ----
    let mut out: Box<dyn Write> = if out_path == "-" {
        Box::new(std::io::stdout())
    } else {
        Box::new(File::create(&out_path)?)
    };

    let cat_json = |name: &str, c: &Cat| -> String {
        let pct_aligned = if c.total_reads > 0 { c.pf_reads_aligned as f64 / c.total_reads as f64 } else { 0.0 };
        let mean_rl = if c.total_reads > 0 { c.read_length_sum as f64 / c.total_reads as f64 } else { 0.0 };
        let pct_in_pairs = if c.pf_reads_aligned > 0 { c.reads_aligned_in_pairs as f64 / c.pf_reads_aligned as f64 } else { 0.0 };
        let strand_balance = if c.pf_reads_aligned > 0 { (c.pf_reads_aligned - c.neg_strand_aligned) as f64 / c.pf_reads_aligned as f64 } else { 0.0 };
        let mismatch_rate = if c.aligned_bases > 0 { c.mismatch_sum as f64 / c.aligned_bases as f64 } else { 0.0 };
        let indel_rate = if c.aligned_bases > 0 { c.indel_events as f64 / c.aligned_bases as f64 } else { 0.0 };
        let hq_error_rate = if c.hq_aligned_bases > 0 { c.hq_mismatch_sum as f64 / c.hq_aligned_bases as f64 } else { 0.0 };
        let chim_den = if c.chim_fragment_den > 0 { c.chim_fragment_den } else { c.reads_aligned_in_pairs };
        let pct_chimeras = if chim_den > 0 { c.chimeric_reads as f64 / chim_den as f64 } else { 0.0 };
        let chim_frac = |x: u64| if chim_den > 0 { x as f64 / chim_den as f64 } else { 0.0 };
        format!(
            "    \"{name}\": {{\n\
             \x20     \"TOTAL_READS\": {},\n\
             \x20     \"PF_READS_ALIGNED\": {},\n\
             \x20     \"PF_ALIGNED_BASES\": {},\n\
             \x20     \"PF_HQ_ALIGNED_READS\": {},\n\
             \x20     \"PF_HQ_ALIGNED_BASES\": {},\n\
             \x20     \"PF_HQ_ALIGNED_Q20_BASES\": {},\n\
             \x20     \"PCT_PF_READS_ALIGNED\": {},\n\
             \x20     \"MEAN_READ_LENGTH\": {},\n\
             \x20     \"READS_ALIGNED_IN_PAIRS\": {},\n\
             \x20     \"PCT_READS_ALIGNED_IN_PAIRS\": {},\n\
             \x20     \"STRAND_BALANCE\": {},\n\
             \x20     \"PF_MISMATCH_RATE\": {},\n\
             \x20     \"PF_HQ_ERROR_RATE\": {},\n\
             \x20     \"PF_INDEL_RATE\": {},\n\
             \x20     \"PCT_CHIMERAS\": {},\n\
             \x20     \"CHIMERAS_BREAKDOWN\": {{\n\
             \x20       \"different_contig\": {}, \"different_contig_pct\": {},\n\
             \x20       \"bad_orientation\": {}, \"bad_orientation_pct\": {},\n\
             \x20       \"large_insert_gt_100kb\": {}, \"large_insert_gt_100kb_pct\": {},\n\
             \x20       \"split_read_sa_tag\": {}, \"split_read_sa_tag_pct\": {},\n\
             \x20       \"bad_orientation_breakdown\": {{\n\
             \x20         \"same_strand_ff\": {}, \"same_strand_ff_pct\": {},\n\
             \x20         \"same_strand_rr\": {}, \"same_strand_rr_pct\": {},\n\
             \x20         \"everted_rf\": {}, \"everted_rf_pct\": {},\n\
             \x20         \"dist_lt_1kb\": {}, \"dist_lt_1kb_pct\": {},\n\
             \x20         \"dist_1_10kb\": {}, \"dist_1_10kb_pct\": {},\n\
             \x20         \"dist_10_100kb\": {}, \"dist_10_100kb_pct\": {},\n\
             \x20         \"dist_gt_100kb\": {}, \"dist_gt_100kb_pct\": {}\n\
             \x20       }}\n\
             \x20     }}\n\
             \x20   }}",
            c.total_reads, c.pf_reads_aligned, c.aligned_bases,
            c.hq_reads, c.hq_aligned_bases, c.hq_q20_bases,
            f6(pct_aligned), f6(mean_rl),
            c.reads_aligned_in_pairs, f6(pct_in_pairs), f6(strand_balance),
            f6(mismatch_rate), f6(hq_error_rate), f6(indel_rate), f6(pct_chimeras),
            c.chim_diff_contig, f6(chim_frac(c.chim_diff_contig)),
            c.chim_bad_orient, f6(chim_frac(c.chim_bad_orient)),
            c.chim_large_insert, f6(chim_frac(c.chim_large_insert)),
            c.chim_split_sa, f6(chim_frac(c.chim_split_sa)),
            c.chim_bo_ff, f6(chim_frac(c.chim_bo_ff)),
            c.chim_bo_rr, f6(chim_frac(c.chim_bo_rr)),
            c.chim_bo_rf, f6(chim_frac(c.chim_bo_rf)),
            c.chim_bo_dist_lt1kb, f6(chim_frac(c.chim_bo_dist_lt1kb)),
            c.chim_bo_dist_1_10kb, f6(chim_frac(c.chim_bo_dist_1_10kb)),
            c.chim_bo_dist_10_100kb, f6(chim_frac(c.chim_bo_dist_10_100kb)),
            c.chim_bo_dist_gt100kb, f6(chim_frac(c.chim_bo_dist_gt100kb)),
        )
    };

    writeln!(out, "{{")?;
    writeln!(out, "  \"sample\": \"{sample}\",")?;
    writeln!(out, "  \"alignment_stat\": {{")?;
    writeln!(out, "{},", cat_json("FIRST_OF_PAIR", &first))?;
    writeln!(out, "{},", cat_json("SECOND_OF_PAIR", &second))?;
    // UNPAIRED (Picard CATEGORY=UNPAIRED) is emitted only when single-end reads were seen,
    // so paired-end output is byte-identical to earlier releases.
    if unpaired.total_reads > 0 {
        writeln!(out, "{},", cat_json("UNPAIRED", &unpaired))?;
    }
    writeln!(out, "{}", cat_json("PAIR", &pair))?;
    writeln!(out, "  }},")?;
    // chimera feature cube (PAIR category): modular substrate for chimera composition.
    // Named mechanism classes are declarative sums over these cells (chimera_classes.yaml).
    {
        let loc = ["same_contig", "diff_contig"];
        let ori = ["proper_FR", "FF", "RR", "everted"];
        let dst = ["na", "lt_1kb", "1_10kb", "10_100kb", "100kb_1Mb", "gt_1Mb"];
        let sav = ["no_SA", "has_SA"];
        writeln!(out, "  \"chimera_cube\": {{")?;
        writeln!(out, "    \"axes\": {{\"locus\": [\"same_contig\", \"diff_contig\"], \"orient\": [\"proper_FR\", \"FF\", \"RR\", \"everted\"], \"dist\": [\"na\", \"lt_1kb\", \"1_10kb\", \"10_100kb\", \"100kb_1Mb\", \"gt_1Mb\"], \"sa\": [\"no_SA\", \"has_SA\"]}},")?;
        writeln!(out, "    \"reads_aligned_in_pairs\": {},", pair.reads_aligned_in_pairs)?;
        write!(out, "    \"cells\": [")?;
        let mut first_cell = true;
        for l in 0..CUBE_LOCUS {
            for o in 0..CUBE_ORIENT {
                for d in 0..CUBE_DIST {
                    for s in 0..CUBE_SA {
                        let cnt = pair.chim_cube.0[cube_idx(l, o, d, s)];
                        if cnt == 0 { continue; }
                        if !first_cell { write!(out, ", ")?; }
                        first_cell = false;
                        write!(out, "{{\"locus\": \"{}\", \"orient\": \"{}\", \"dist\": \"{}\", \"sa\": \"{}\", \"count\": {}}}",
                               loc[l], ori[o], dst[d], sav[s], cnt)?;
                    }
                }
            }
        }
        writeln!(out, "]")?;
        writeln!(out, "  }},")?;
    }
    // chrM fraction denominator: all aligned bases (paired + single-end). chrm_aligned_bases
    // already includes single-end reads, so the denominator must too. Unchanged for
    // paired-only input (unpaired.aligned_bases == 0).
    let total_aln = pair.aligned_bases + unpaired.aligned_bases;
    let pct_chrm = if total_aln > 0 { chrm_aligned_bases as f64 / total_aln as f64 } else { 0.0 };
    writeln!(out, "  \"contig\": {{")?;
    writeln!(out, "    \"cov_total_bases\": {},", total_aln)?;
    writeln!(out, "    \"cov_chrm_bases\": {},", chrm_aligned_bases)?;
    writeln!(out, "    \"pct_chrm\": {}", f6(pct_chrm))?;
    writeln!(out, "  }},")?;
    writeln!(out, "  \"optical\": {{")?;
    writeln!(out, "    \"read_pair_optical_duplicates\": {}", optical_dups)?;
    writeln!(out, "  }},")?;
    writeln!(out, "  \"preseq\": {{")?;
    match preseq_count {
        Some(v) => writeln!(out, "    \"preseq_count\": {},", v.round() as i64)?,
        None => writeln!(out, "    \"preseq_count\": null,")?,
    }
    writeln!(out, "    \"total_bins\": {},", preseq_total_bins as i64)?;
    writeln!(out, "    \"distinct_bins\": {},", preseq_distinct_bins as i64)?;
    writeln!(out, "    \"max_coverage\": {},", preseq_max_cov)?;
    writeln!(out, "    \"counts_of_1\": {}", *preseq_acc.hist.get(1).unwrap_or(&0.0) as i64)?;
    writeln!(out, "  }},")?;
    // Lorenz-curve / Gini coefficient (bam-lorenz-coverage reproduction).
    writeln!(out, "  \"lorenz\": {{")?;
    match lorenz_roc {
        Some((roc, cov_pos, bases)) => {
            writeln!(out, "    \"gini_coefficient_index\": {},", f9(0.5 - roc))?;
            writeln!(out, "    \"roc_lorenz_curve\": {},", f9(roc))?;
            writeln!(out, "    \"total_covered_positions_of_genome\": {},", cov_pos)?;
            writeln!(out, "    \"total_sequenced_bases\": {},", bases)?;
            writeln!(out, "    \"total_investigated_genomic_positions\": {}", lorenz.investigated)?;
        }
        None => {
            writeln!(out, "    \"gini_coefficient_index\": null,")?;
            writeln!(out, "    \"roc_lorenz_curve\": null,")?;
            writeln!(out, "    \"total_covered_positions_of_genome\": 0,")?;
            writeln!(out, "    \"total_sequenced_bases\": 0,")?;
            writeln!(out, "    \"total_investigated_genomic_positions\": {}", lorenz.investigated)?;
        }
    }
    writeln!(out, "  }},")?;
    let qy_q20_rate = if qy_total_bases > 0 { qy_q20_bases as f64 / qy_total_bases as f64 } else { 0.0 };
    let qy_q30_rate = if qy_total_bases > 0 { qy_q30_bases as f64 / qy_total_bases as f64 } else { 0.0 };
    writeln!(out, "  \"quality_yield\": {{")?;
    writeln!(out, "    \"TOTAL_READS\": {},", qy_reads)?;
    writeln!(out, "    \"TOTAL_BASES\": {},", qy_total_bases)?;
    writeln!(out, "    \"Q20_BASES\": {},", qy_q20_bases)?;
    writeln!(out, "    \"Q30_BASES\": {},", qy_q30_bases)?;
    writeln!(out, "    \"PCT_Q20\": {},", f6(qy_q20_rate))?;
    writeln!(out, "    \"PCT_Q30\": {},", f6(qy_q30_rate))?;
    // QUALITY_HISTOGRAM[q] = bases observed at phred q, over exactly the bases counted in
    // TOTAL_BASES. Consumed by the theoretical HET_SNP_SENSITIVITY model.
    writeln!(out, "    \"QUALITY_HISTOGRAM\": [{}]",
        qual_hist.iter().map(|c| c.to_string()).collect::<Vec<_>>().join(","))?;
    writeln!(out, "  }},")?;
    writeln!(out, "  \"insert_size\": {{")?;
    writeln!(out, "    \"MEDIAN_INSERT_SIZE\": {},", f6(ins_median))?;
    writeln!(out, "    \"MEDIAN_ABSOLUTE_DEVIATION\": {},", f6(ins_mad))?;
    writeln!(out, "    \"MIN_INSERT_SIZE\": {},", ins_min)?;
    writeln!(out, "    \"MAX_INSERT_SIZE\": {},", ins_max)?;
    writeln!(out, "    \"MEAN_INSERT_SIZE\": {},", f6(ins_mean))?;
    writeln!(out, "    \"STANDARD_DEVIATION\": {},", f6(ins_sd))?;
    writeln!(out, "    \"READ_PAIRS\": {}", ins_pairs)?;
    if cov.enabled {
        writeln!(out, "  }},")?;
        let b = (cov.counted_bases + cov.exc_mapq + cov.exc_dupe + cov.exc_unpaired
            + cov.exc_baseq + cov.exc_overlap) as f64; // total aligned bases; capped bases already counted
        let excf = |x: u64| if b > 0.0 { x as f64 / b } else { 0.0 };
        let mean = cov.mean();
        let exc_total = cov.exc_mapq + cov.exc_dupe + cov.exc_unpaired + cov.exc_baseq
            + cov.exc_overlap + cov.exc_capped;
        writeln!(out, "  \"coverage\": {{")?;
        writeln!(out, "    \"GENOME_TERRITORY\": {},", cov.territory)?;
        let cov_median = cov.median();
        writeln!(out, "    \"MEAN_COVERAGE\": {},", f6(mean))?;
        writeln!(out, "    \"SD_COVERAGE\": {},", f6(cov.sd(mean)))?;
        writeln!(out, "    \"MEDIAN_COVERAGE\": {},", f6(cov_median))?;
        writeln!(out, "    \"MAD_COVERAGE\": {},", f6(cov.mad(cov_median)))?;
        writeln!(out, "    \"PCT_EXC_MAPQ\": {},", f6(excf(cov.exc_mapq)))?;
        writeln!(out, "    \"PCT_EXC_DUPE\": {},", f6(excf(cov.exc_dupe)))?;
        writeln!(out, "    \"PCT_EXC_UNPAIRED\": {},", f6(excf(cov.exc_unpaired)))?;
        writeln!(out, "    \"PCT_EXC_BASEQ\": {},", f6(excf(cov.exc_baseq)))?;
        writeln!(out, "    \"PCT_EXC_OVERLAP\": {},", f6(excf(cov.exc_overlap)))?;
        writeln!(out, "    \"PCT_EXC_CAPPED\": {},", f6(excf(cov.exc_capped)))?;
        writeln!(out, "    \"PCT_EXC_TOTAL\": {},", f6(excf(exc_total)))?;
        writeln!(out, "    \"PCT_1X\": {},", f6(cov.pct_ge(1)))?;
        writeln!(out, "    \"PCT_5X\": {},", f6(cov.pct_ge(5)))?;
        writeln!(out, "    \"PCT_10X\": {},", f6(cov.pct_ge(10)))?;
        writeln!(out, "    \"PCT_15X\": {},", f6(cov.pct_ge(15)))?;
        writeln!(out, "    \"PCT_20X\": {},", f6(cov.pct_ge(20)))?;
        writeln!(out, "    \"PCT_25X\": {},", f6(cov.pct_ge(25)))?;
        writeln!(out, "    \"PCT_30X\": {},", f6(cov.pct_ge(30)))?;
        writeln!(out, "    \"PCT_40X\": {},", f6(cov.pct_ge(40)))?;
        writeln!(out, "    \"PCT_50X\": {},", f6(cov.pct_ge(50)))?;
        writeln!(out, "    \"PCT_60X\": {},", f6(cov.pct_ge(60)))?;
        writeln!(out, "    \"PCT_70X\": {},", f6(cov.pct_ge(70)))?;
        writeln!(out, "    \"PCT_80X\": {},", f6(cov.pct_ge(80)))?;
        writeln!(out, "    \"PCT_90X\": {},", f6(cov.pct_ge(90)))?;
        writeln!(out, "    \"PCT_100X\": {},", f6(cov.pct_ge(100)))?;
        // DEPTH_HISTOGRAM[d] = genomic positions at depth d, d = 0..=coverage_cap, summing
        // to GENOME_TERRITORY. cov.hist only holds covered positions (finalize_pos is never
        // called with depth 0), so index 0 is reconstructed as the uncovered remainder --
        // the same identity median() uses. Feeds the HET_SNP_SENSITIVITY depth distribution.
        let uncovered = cov.territory.saturating_sub(cov.covered_positions);
        writeln!(out, "    \"DEPTH_HISTOGRAM\": [{}]",
            std::iter::once(uncovered)
                .chain(cov.hist.iter().skip(1).copied())
                .map(|c| c.to_string())
                .collect::<Vec<_>>()
                .join(","))?;
        writeln!(out, "  }}{}", if gcbias.is_some() { "," } else { "" })?;
    } else {
        writeln!(out, "  }}{}", if gcbias.is_some() { "," } else { "" })?;
    }
    if let Some(ref gb) = gcbias {
        let total_windows: u64 = gb.windows_by_gc.iter().sum();
        let total_reads: u64 = gb.reads_by_gc.iter().sum();
        let mean_rpw = if total_windows > 0 { total_reads as f64 / total_windows as f64 } else { 0.0 };
        let mut at_dropout = 0.0f64;
        let mut gc_dropout = 0.0f64;
        for i in 0..101usize {
            let rr = if total_reads > 0 { gb.reads_by_gc[i] as f64 / total_reads as f64 } else { 0.0 };
            let rw = if total_windows > 0 { gb.windows_by_gc[i] as f64 / total_windows as f64 } else { 0.0 };
            let d = (rw - rr) * 100.0;
            if d > 0.0 {
                if i <= 50 { at_dropout += d; } else { gc_dropout += d; }
            }
        }
        writeln!(out, "  \"gcbias\": {{")?;
        writeln!(out, "    \"WINDOW_SIZE\": {},", gb.window)?;
        writeln!(out, "    \"TOTAL_CLUSTERS\": {},", gb.total_clusters)?;
        writeln!(out, "    \"ALIGNED_READS\": {},", gb.total_aligned_reads)?;
        writeln!(out, "    \"AT_DROPOUT\": {},", f6(at_dropout))?;
        writeln!(out, "    \"GC_DROPOUT\": {},", f6(gc_dropout))?;
        writeln!(out, "    \"GC_NC_0_19\": {},", f6(gb.norm_cov(mean_rpw, 0, 19)))?;
        writeln!(out, "    \"GC_NC_20_39\": {},", f6(gb.norm_cov(mean_rpw, 20, 39)))?;
        writeln!(out, "    \"GC_NC_40_59\": {},", f6(gb.norm_cov(mean_rpw, 40, 59)))?;
        writeln!(out, "    \"GC_NC_60_79\": {},", f6(gb.norm_cov(mean_rpw, 60, 79)))?;
        writeln!(out, "    \"GC_NC_80_100\": {}", f6(gb.norm_cov(mean_rpw, 80, 100)))?;
        writeln!(out, "  }}")?;
    }
    writeln!(out, "}}")?;

    Ok(())
}

/// Dense insert-size histogram limit. Larger TLEN values are retained exactly in a
/// sparse overflow map so median/MAD and extrema remain canonical without unbounded RAM.
const INSERT_CAP: usize = 1_000_000;

/// Base-quality histogram length: phred 0..=93 inclusive (the SAM spec maximum, ASCII '~'
/// minus 33). Feeds the theoretical HET_SNP_SENSITIVITY model, which needs the full
/// quality distribution rather than the Q20/Q30 threshold counts.
const QUAL_HIST_LEN: usize = 94;

fn insert_values<'a>(
    hist: &'a [u64],
    overflow: &'a std::collections::BTreeMap<u64, u64>,
) -> impl Iterator<Item = (u64, u64)> + 'a {
    hist.iter()
        .enumerate()
        .filter_map(|(value, &count)| (count > 0).then_some((value as u64, count)))
        .chain(overflow.iter().map(|(&value, &count)| (value, count)))
}

/// 0-indexed k-th value from the exact insert-size distribution.
fn hist_kth(
    hist: &[u64],
    overflow: &std::collections::BTreeMap<u64, u64>,
    k: u64,
) -> u64 {
    let mut cumulative = 0u64;
    for (value, count) in insert_values(hist, overflow) {
        cumulative += count;
        if cumulative > k {
            return value;
        }
    }
    overflow
        .keys()
        .next_back()
        .copied()
        .unwrap_or_else(|| hist.len().saturating_sub(1) as u64)
}

/// Insert-size stats matching Picard CollectInsertSizeMetrics: MEDIAN over all pairs,
/// then MEAN/SD computed only over pairs within MEDIAN +/- DEVIATIONS*MAD (default 10).
fn insert_stats_hist(
    hist: &[u64],
    overflow: &std::collections::BTreeMap<u64, u64>,
) -> (f64, f64, f64, f64, u64) {
    let n: u64 = insert_values(hist, overflow).map(|(_, count)| count).sum();
    if n == 0 {
        return (0.0, 0.0, 0.0, 0.0, 0);
    }

    let median = if n % 2 == 1 {
        hist_kth(hist, overflow, n / 2) as f64
    } else {
        (hist_kth(hist, overflow, n / 2 - 1) as f64
            + hist_kth(hist, overflow, n / 2) as f64)
            / 2.0
    };

    // Store absolute deviations in half-base units so even-count medians retain .5.
    let mut deviation_counts: std::collections::BTreeMap<u64, u64> =
        std::collections::BTreeMap::new();
    for (value, count) in insert_values(hist, overflow) {
        let deviation_twice = ((value as f64 - median).abs() * 2.0).round() as u64;
        *deviation_counts.entry(deviation_twice).or_insert(0) += count;
    }
    let deviation_kth = |k: u64| -> u64 {
        let mut cumulative = 0u64;
        for (&deviation_twice, &count) in deviation_counts.iter() {
            cumulative += count;
            if cumulative > k {
                return deviation_twice;
            }
        }
        deviation_counts.keys().next_back().copied().unwrap_or(0)
    };
    let mad = if n % 2 == 1 {
        deviation_kth(n / 2) as f64 / 2.0
    } else {
        (deviation_kth(n / 2 - 1) as f64 + deviation_kth(n / 2) as f64) / 4.0
    };

    const DEVIATIONS: f64 = 10.0;
    let lower = median - DEVIATIONS * mad;
    let upper = median + DEVIATIONS * mad;
    let mut retained = 0u64;
    let mut sum = 0.0;
    for (value, count) in insert_values(hist, overflow) {
        let x = value as f64;
        if x >= lower && x <= upper {
            retained += count;
            sum += x * count as f64;
        }
    }
    if retained == 0 {
        return (median, mad, median, 0.0, n);
    }

    let mean = sum / retained as f64;
    let mut squared_deviation_sum = 0.0;
    for (value, count) in insert_values(hist, overflow) {
        let x = value as f64;
        if x >= lower && x <= upper {
            let deviation = x - mean;
            squared_deviation_sum += deviation * deviation * count as f64;
        }
    }
    let variance = squared_deviation_sum / retained as f64;
    (median, mad, mean, variance.sqrt(), n)
}

// Extract an unsigned integer from a BAM data field value (NM is typically u8/u16/i32).
fn value_as_u64(
    value: &noodles_sam::alignment::record::data::field::Value,
) -> Option<u64> {
    use noodles_sam::alignment::record::data::field::Value;
    match value {
        Value::Int8(v) => Some(*v as u64),
        Value::UInt8(v) => Some(*v as u64),
        Value::Int16(v) => Some(*v as u64),
        Value::UInt16(v) => Some(*v as u64),
        Value::Int32(v) => Some(*v as u64),
        Value::UInt32(v) => Some(*v as u64),
        _ => None,
    }
}

fn value_as_string(
    value: &noodles_sam::alignment::record::data::field::Value,
) -> Option<String> {
    use noodles_sam::alignment::record::data::field::Value;
    match value {
        Value::String(s) => Some(String::from_utf8_lossy(s.as_ref()).into_owned()),
        _ => None,
    }
}

/// Fallback substitution estimate when MD is absent: NM edit distance minus indel bases.
fn nm_minus_indels(record: &bam::Record, _aln_bases: u64) -> u64 {
    use noodles_sam::alignment::record::cigar::op::Kind;
    let mut indel_bases = 0u64;
    for op in record.cigar().iter().flatten() {
        if matches!(op.kind(), Kind::Insertion | Kind::Deletion) {
            indel_bases += op.len() as u64;
        }
    }
    if let Some(Ok(v)) = record.data().get(&Tag::from(*b"NM")) {
        if let Some(nm) = value_as_u64(&v) {
            return nm.saturating_sub(indel_bases);
        }
    }
    0
}
