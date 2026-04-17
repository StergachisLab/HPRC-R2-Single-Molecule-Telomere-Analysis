#!/usr/bin/env python3
# mp_telomeric_meanq_filter_bam.py
# Parallel telomeric mean-Q filter that WRITES BAMs (per-chrom, then merged).
#
# Usage:
#   python3 mp_telomeric_meanq_filter_bam.py reads.bam telomeres.bed \
#       --min-meanq 20 -p 8 -o telomeric_q20.bam --merge-threads 8
#
import argparse, os, sys, tempfile, multiprocessing as mp, subprocess, shutil
from bisect import bisect_left
import pysam

# CIGAR op codes
M,I,D,N,S,H,P,EQ,X = 0,1,2,3,4,5,6,7,8

def merge_intervals(iv):
    if not iv: return []
    iv = sorted(iv)
    out = [iv[0]]
    for s,e in iv[1:]:
        ls,le = out[-1]
        if s <= le:
            if e > le: out[-1] = (ls, e)
        else:
            out.append((s,e))
    return out

def load_telomeres_bed(path):
    # BED: chrom start end [name] [score] [strand]; we ignore strand for filtering
    chrom2ivs = {}
    with open(path) as f:
        for ln in f:
            if not ln.strip() or ln.startswith("#"): continue
            parts = ln.rstrip().split("\t")
            if len(parts) < 3: continue
            chrom = parts[0]; s = int(parts[1]); e = int(parts[2])
            chrom2ivs.setdefault(chrom, []).append((s,e))
    chrom2ivs = {c: merge_intervals(v) for c,v in chrom2ivs.items()}
    chrom2starts = {c: [s for s,_ in ivs] for c,ivs in chrom2ivs.items()}
    return chrom2ivs, chrom2starts

def ref_in_iv(chrom, pos, ivs, starts):
    arr = starts.get(chrom)
    if not arr: return -1
    i = bisect_left(arr, pos+1) - 1
    if i < 0: return -1
    s,e = ivs[chrom][i]
    return i if (s <= pos < e) else -1

def sumq_span(q, qpos, length):
    # fast enough on Python; q is an array of ints
    return sum(q[qpos:qpos+length])

def telomere_meanQ_for_read(aln, ivs, starts):
    """
    Compute mean PHRED by telomere interval(s) overlapped by read bases.
    Include: M,=,X; I; end S if adjacent ref edge is telomeric. Exclude D,N.
    Returns (means:list[float], counted:bool)
    """
    q = aln.query_qualities
    if q is None: return [], False
    chrom = aln.reference_name
    if chrom not in ivs: return [], False
    ivlist = ivs[chrom]

    sumq = {}
    cnt  = {}

    rpos = aln.reference_start
    qpos = 0
    ct   = aln.cigartuples or []

    # Left soft-clip: include if first aligned ref base is telomeric
    if ct and ct[0][0] == S:
        lS = ct[0][1]
        if lS > 0 and rpos is not None:
            bi = ref_in_iv(chrom, rpos, ivs, starts)
            if bi >= 0:
                sumq[bi] = sumq.get(bi,0) + sumq_span(q, 0, lS)
                cnt[bi]  = cnt.get(bi,0)  + lS
        qpos += lS

    # Core ops
    for op,len_ in ct:
        if op in (H, S):  # end S handled above/below
            continue
        if op in (M, EQ, X):
            seg_start = rpos
            seg_end   = rpos + len_
            # start scanning intervals overlapping this segment
            idx = bisect_left(starts[chrom], seg_start+1) - 1
            if idx < 0: idx = 0
            while idx < len(ivlist):
                s,e = ivlist[idx]
                if s >= seg_end: break
                if e > seg_start and s < seg_end:
                    lo = max(seg_start, s)
                    hi = min(seg_end,   e)
                    L  = hi - lo
                    offs = (lo - seg_start)
                    sumq[idx] = sumq.get(idx,0) + sumq_span(q, qpos + offs, L)
                    cnt[idx]  = cnt.get(idx,0)  + L
                idx += 1
            rpos += len_
            qpos += len_
        elif op == I:
            if len_ > 0 and rpos is not None:
                bi = ref_in_iv(chrom, rpos, ivs, starts)
                if bi >= 0:
                    sumq[bi] = sumq.get(bi,0) + sumq_span(q, qpos, len_)
                    cnt[bi]  = cnt.get(bi,0)  + len_
                qpos += len_
        elif op in (D, N):
            rpos += len_
        else:
            pass

    # Right soft-clip: include if last aligned ref base (rpos-1) is telomeric
    if ct and ct[-1][0] == S:
        rS = ct[-1][1]
        if rS > 0 and rpos is not None:
            edge = rpos - 1
            bi = ref_in_iv(chrom, edge, ivs, starts)
            if bi >= 0:
                sumq[bi] = sumq.get(bi,0) + sumq_span(q, qpos, rS)
                cnt[bi]  = cnt.get(bi,0)  + rS
        qpos += rS

    means = []
    counted = False
    for bi in sumq:
        if cnt[bi] > 0:
            means.append(sumq[bi] / float(cnt[bi]))
            counted = True
    return means, counted

def worker(args):
    chrom, bam_path, telos_iv, telos_st, min_meanq, min_mapq, flag_filter, all_sides, out_bam_path = args
    kept = tot = 0
    with pysam.AlignmentFile(bam_path, "rb") as bam, \
         pysam.AlignmentFile(out_bam_path, "wb", header=bam.header) as outbam:
        if chrom not in telos_iv:
            return out_bam_path, kept, tot
        for aln in bam.fetch(chrom):
            tot += 1
            if aln.flag & flag_filter: continue
            if aln.mapping_quality < min_mapq: continue
            means, counted = telomere_meanQ_for_read(aln, telos_iv, telos_st)
            if not counted: continue
            keep = all(m >= min_meanq for m in means) if all_sides else any(m >= min_meanq for m in means)
            if keep:
                outbam.write(aln)
                kept += 1
    return out_bam_path, kept, tot

def run_cmd(cmd):
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\nSTDERR:\n{proc.stderr}")
    return proc

def main():
    ap = argparse.ArgumentParser(description="Parallel telomeric mean-Q filter that writes BAM directly (per-chrom + merged).")
    ap.add_argument("bam", help="Input BAM (indexed)")
    ap.add_argument("telomere_bed", help="BED with telomeric intervals (strand +/- optional)")
    ap.add_argument("--min-meanq", type=float, required=True, help="Minimum mean PHRED across telomeric read bases")
    ap.add_argument("--min-mapq", type=int, default=0, help="Mapping quality cutoff (default 0)")
    ap.add_argument("--flag-filter", type=int, default=0x4, help="Exclude flags mask (default drops unmapped)")
    group = ap.add_mutually_exclusive_group()
    group.add_argument("--all-sides", action="store_true", help="Require ALL telomere intervals overlapped by a read to meet threshold")
    group.add_argument("--any-side", action="store_true", help="Keep if ANY telomere interval meets threshold (default)")
    ap.add_argument("-p", "--procs", type=int, default=0, help="Workers (default: min(#chrom_with_telomere, CPU count))")
    ap.add_argument("-o", "--out-bam", default="telomeric_filtered.bam", help="Final merged BAM")
    ap.add_argument("--merge-threads", type=int, default=4, help="Threads for samtools merge")
    ap.add_argument("--tmpdir", default=None, help="Temporary directory (kept for inspection if provided)")
    args = ap.parse_args()
    any_side = (not args.all_sides)

    telos_iv, telos_st = load_telomeres_bed(args.telomere_bed)

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        chroms = [c for c in bam.references if c in telos_iv]

    if not chroms:
        print("No chromosomes in BAM have telomere intervals; nothing to do.", file=sys.stderr)
        # write empty BAM with header copied
        with pysam.AlignmentFile(args.bam, "rb") as ib, pysam.AlignmentFile(args.out_bam, "wb", header=ib.header) as ob:
            pass
        pysam.index(args.out_bam)
        return

    ncpu = os.cpu_count() or 1
    nprocs = args.procs if args.procs and args.procs > 0 else min(len(chroms), ncpu)

    tmpdir = args.tmpdir or tempfile.mkdtemp(prefix="teloq_bam_")

    print(f"[info] processing {len(chroms)} chroms with {nprocs} workers; tmpdir={tmpdir}", file=sys.stderr)

    # plan work: one per chrom, to distinct BAM paths
    work = []
    for c in chroms:
        out_bam_path = os.path.join(tmpdir, f"{c}.keep.bam")
        work.append((c, args.bam, telos_iv, telos_st, args.min_meanq, args.min_mapq,
                     args.flag_filter, (not any_side), out_bam_path))

    kept_total = tot_total = 0
    out_parts = []
    with mp.Pool(processes=nprocs) as pool:
        for out_bam_path, kept, tot in pool.imap_unordered(worker, work):
            out_parts.append(out_bam_path)
            kept_total += kept
            tot_total  += tot
            print(f"[done] {os.path.basename(out_bam_path)}  kept {kept}/{tot}", file=sys.stderr)

    # samtools merge (coordinate-sorted inputs are not required; merge can sort if -n/-t not used; we’ll fill RG/headers from first)
    # For large sets, samtools merge is faster & memory-savvy
    merge_cmd = ["samtools", "merge", "-f", "-@", str(args.merge_threads), args.out_bam] + out_parts
    run_cmd(merge_cmd)
    # coordinate sort to be safe, then index
    sort_prefix = os.path.join(tmpdir, "sorted")
    run_cmd(["samtools", "sort", "-@", str(args.merge_threads), "-o", args.out_bam, args.out_bam])
    run_cmd(["samtools", "index", args.out_bam])

    print(f"[summary] kept {kept_total} / {tot_total} records -> {args.out_bam}", file=sys.stderr)

    if args.tmpdir is None:
        # auto-clean tmpdir
        try:
            shutil.rmtree(tmpdir)
        except Exception:
            pass
    else:
        print(f"[note] temp files kept at: {tmpdir}", file=sys.stderr)

if __name__ == "__main__":
    main()

