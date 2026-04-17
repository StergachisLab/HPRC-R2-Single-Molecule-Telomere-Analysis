#!/usr/bin/env python3
import numpy as np
import pandas as pd
from glob import glob
from tqdm import tqdm
import matplotlib.pyplot as plt
import re
import sys
import warnings
import os
from concurrent.futures import ThreadPoolExecutor, as_completed  # multithreading
from concurrent.futures import ProcessPoolExecutor, as_completed


def get_label_from_series(series: pd.Series) -> str:
    if series.empty or series.isna().all():
        return "unknown"
    first_val = series.dropna().iloc[0]
    parts = str(first_val).split("#")
    if len(parts) >= 2:
        return f"{parts[0]}_{parts[1]}"
    return parts[0]

def rotations(s: str):
    return [s[i:] + s[:i] for i in range(len(s))]


def ID_TVR_PHASE_DISRUPTION(s):
    plus6 = rotations("CCCTAA")
    minus6 = rotations("TTAGGG")
    pat_plus6 = re.compile(r"(?=(?:%s))" % "|".join(plus6), re.I)
    pat_minus6 = re.compile(r"(?=(?:%s))" % "|".join(minus6), re.I)
    plus_idx = np.fromiter((m.start() for m in pat_plus6.finditer(s)), dtype=int)
    minus_idx = np.fromiter((m.start() for m in pat_minus6.finditer(s)), dtype=int)

    is_plus = (plus_idx.size >= minus_idx.size)
    true_idx = plus_idx if is_plus else minus_idx
    if true_idx.size == 0:
        return np.empty((0, 2), dtype=int)

    L = len(s)
    diff = np.zeros(L + 1, dtype=np.int32)
    np.add.at(diff, true_idx, 1)
    np.add.at(diff, true_idx + 6, -1)
    cov = np.cumsum(diff)[:L].astype(np.uint16)

    canon = "CCCTAA" if is_plus else "TTAGGG"
    pat_canon = re.compile(rf"(?={canon})", re.I)
    canon_idx = np.fromiter((m.start() for m in pat_canon.finditer(s)), dtype=int)
    if canon_idx.size:
        diff2 = np.zeros(L + 1, dtype=np.int16)
        np.add.at(diff2, canon_idx, 1)
        np.add.at(diff2, canon_idx + 6, -1)
        canon_mask = np.cumsum(diff2)[:L] > 0
        cov[canon_mask] = 6

    low = cov < 6
    d = np.diff(np.pad(low.astype(np.int8), (1, 1)))
    starts = np.where(d == 1)[0]
    ends = np.where(d == -1)[0]
    spans = np.column_stack([starts, ends - starts]).astype(int)
    return spans

def assignTVRAndAddTelomereSequenceColumn(tsv,
                                          minimum_cutoff_for_read_length=2_500,
                                          maximum_cutoff_for_tvr_length=300):
    teloseqs, tvrs = [], []
    for qn, grp in tsv.groupby("query_name", sort=False):
        seq = grp["query_sequence"].iloc[0]
        qe = int(grp["centered_query_end"].iloc[0])
        if not isinstance(seq, str) or qe < minimum_cutoff_for_read_length:
            teloseq = ""
            tvr = np.empty((0, 2), int)
        else:
            qe = min(qe, len(seq))
            teloseq = seq[-qe:]
            tvr = ID_TVR_PHASE_DISRUPTION(teloseq)
            tvr = tvr[tvr[:, 1] <= maximum_cutoff_for_tvr_length] if tvr.size else tvr
        teloseqs.extend([teloseq] * len(grp))
        tvrs.extend([tvr] * len(grp))
    tsv = tsv.assign(teloseq=teloseqs, tvr=tvrs)
    return tsv


def consensus_telomere_from_group(
    group,
    seq_col: str = "query_sequence",
    end_len_col: str = "centered_query_end",
) -> str:
    tails = []
    for seq, qend in zip(group[seq_col].values, group[end_len_col].values):
        if isinstance(seq, str):
            try:
                q = int(qend)
            except Exception:
                continue
            if 0 < q <= len(seq):
                tails.append(seq[-q:])
    if not tails:
        return ""
    Lmax = max(len(s) for s in tails)
    medianL = int(np.nanmedian([len(s) for s in tails]))
    arr = np.full((len(tails), Lmax), '-', dtype='<U1')
    for i, s in enumerate(tails):
        arr[i, :len(s)] = list(s)
    bases = np.array(list("ACGT"))
    counts = np.stack([(arr == b).sum(axis=0) for b in bases], axis=0)
    winners = counts.argmax(axis=0)
    maxcounts = counts.max(axis=0)
    cons = bases[winners]
    cons = np.where(maxcounts == 0, '-', cons)
    return ''.join(cons)[:medianL]

def collect_consensus_telomeres(
    all_tsvs,
    chrom_end: str,
    min_reads_per_group: int = 3,
    only_labels=None,
    t_cap=2000
):
    if isinstance(only_labels, str):
        only_labels = {only_labels}
    elif only_labels is not None:
        only_labels = set(only_labels)
    rows = []
    for tsv in all_tsvs:
        tmp = tsv[tsv.chrom_assignment == chrom_end]
        if tmp.empty:
            continue
        grouped = (
            tmp.assign(span=tmp["centered_query_end"].astype(int))
               .sort_values("span", ascending=True)
               .groupby(["chrom", "centering_position"])
        )
        for (chrom, center_pos), group in grouped:
            if len(group) < min_reads_per_group:
                continue
            label = get_label_from_series(group["chrom"])
            if (only_labels is not None) and (label not in only_labels):
                continue
            cons = consensus_telomere_from_group(group, "query_sequence", "centered_query_end")
            if not isinstance(cons, str) or len(cons) == 0:
                continue
            if t_cap is not None:
                cons = cons[:t_cap]
            rows.append({
                "label": label,
                "assignment": chrom_end,
                "contig": chrom,
                "centering_position": int(center_pos),
                "consensus": cons,
                "length_bp": len(cons),
                "n_reads": int(len(group)),
            })
    return pd.DataFrame(rows).sort_values(
        ["label", "assignment", "centering_position"]
    ).reset_index(drop=True)


def _encode_seq_2bit(seq: str) -> np.ndarray:
    seq = (seq or "").upper()
    lut = np.full(128, 255, dtype=np.uint8)
    lut[ord('A')] = 0; lut[ord('C')] = 1; lut[ord('G')] = 2; lut[ord('T')] = 3
    b = np.frombuffer(seq.encode('ascii', 'ignore'), dtype=np.uint8)
    return lut[b]

def _roll_kmer_codes(xs: np.ndarray, k: int) -> np.ndarray:
    if k <= 0 or xs.size == 0:
        return np.zeros(0, dtype=np.uint64)
    mask = (1 << (2 * k)) - 1
    codes = []
    code = 0
    run = 0
    for v in xs:
        if v < 4:  # valid A/C/G/T
            code = ((code << 2) | int(v)) & mask
            run += 1
            if run >= k:
                codes.append(code)
        else:
            code = 0
            run = 0
    return np.array(codes, dtype=np.uint64)

def _canon_codes(k: int):
    canon = ["TTAGGG", "CCCTAA"]
    codes = set()
    for km in canon:
        xs = _encode_seq_2bit(km[:k].upper())
        code = 0
        for b in xs:
            code = (code << 2) | int(b)
        codes.add(code)
    return codes

_CANON_CACHE = {}

def _get_canon_set(k):
    cs = _CANON_CACHE.get(k)
    if cs is None:
        cs = _canon_codes(k)
        _CANON_CACHE[k] = cs
    return cs

def sw_row_update(prev_row, bk, bk_is_canon, ai, ai_is_canon, match, mismatch, gap, tvrs_only):

    nb = bk.size
    s_row = np.where(bk == ai, match, mismatch)
    if tvrs_only and (ai_is_canon or bk_is_canon is not None):
        s_row = np.where((bk_is_canon) | ai_is_canon, 0.0, s_row)

    diag = prev_row[:-1] + s_row
    up   = prev_row[1:] + gap
    temp = np.maximum(diag, up)

    g = gap
    z = temp - np.arange(1, nb + 1) * g
    left = np.maximum.accumulate(z) + np.arange(1, nb + 1) * g

    cur = np.maximum(0.0, np.maximum(temp, left))
    cur_row = np.empty(nb + 1, dtype=float)
    cur_row[0] = 0.0
    cur_row[1:] = cur
    return cur_row

def kmer_sw_similarity_fast(a: str, b: str, k: int = 6, match: float = 2.0,
                            mismatch: float = -1.0, gap: float = -2.0,
                            tvrs_only: bool = True, return_dp: bool = False):
    A = _encode_seq_2bit(a); B = _encode_seq_2bit(b)
    ak = _roll_kmer_codes(A, k); bk = _roll_kmer_codes(B, k)
    na, nb = ak.size, bk.size
    if na == 0 or nb == 0:
        return (np.zeros((1,1)), na, nb) if return_dp else 0.0

    canon = _get_canon_set(k) if tvrs_only else None
    ak_is_canon = np.fromiter((c in canon for c in ak), dtype=bool, count=na) if canon else None
    bk_is_canon = np.fromiter((c in canon for c in bk), dtype=bool, count=nb) if canon else None

    H = np.zeros((na + 1, nb + 1), dtype=float)
    P = np.zeros_like(H)  # prefix-max matrix
    best = 0.0
    for i in range(1, na + 1):
        ai = ak[i-1]
        ai_c = ak_is_canon[i-1] if canon is not None else False
        cur = sw_row_update(H[i-1], bk, bk_is_canon, ai, ai_c, match, mismatch, gap, tvrs_only)
        H[i] = cur
        P[i, 1:] = np.maximum(P[i-1, 1:], np.maximum.accumulate(H[i, 1:]))
        best = max(best, P[i, nb])

    if return_dp:
        return P, na, nb
    ideal = max(1.0, match * min(na, nb))
    return min(1.0, best / ideal)

def kmer_sw_sim_for_tcap(P: np.ndarray, na: int, nb: int, t_cap_a: int, t_cap_b: int,
                         k: int = 6, match: float = 2.0) -> float:
    ia = max(0, min(na, t_cap_a - k + 1))
    jb = max(0, min(nb, t_cap_b - k + 1))
    if ia == 0 or jb == 0:
        return 0.0
    best = P[ia, jb]
    ideal = max(1.0, match * min(ia, jb))
    return min(1.0, best / ideal)


def pairwise_tvr_similarity_alignment(
    df_consensus, t_cap_range,
    min_tvr_size: int = 1,
    max_tvr_size: int = 300,
    k: int = 6,
):
    out = []
    for assignment, g in df_consensus.groupby('assignment'):
        rows = []
        for _, row in g.iterrows():
            seq = row["consensus"]
            rows.append((row["label"], int(row["length_bp"]), seq, row["contig"]))
        for i in range(len(rows)):
            li, ti, si, contig_i = rows[i]
            for j in range(i + 1, len(rows)):
                lj, tj, sj, contig_j = rows[j]
                if not si or not sj:
                    sims_by_cap = {}
                else:
                    P, na, nb = kmer_sw_similarity_fast(si, sj, k=k, tvrs_only=True, return_dp=True)
                    sims_by_cap = {t_cap: kmer_sw_sim_for_tcap(P, na, nb, t_cap, t_cap, k=k)
                                   for t_cap in t_cap_range}
                out.append({
                    "contig_i": contig_i,
                    "contig_j": contig_j,
                    "assignment": assignment,
                    "label_i": li,
                    "label_j": lj,
                    "sim": sims_by_cap,
                })
    df = pd.DataFrame(out).sort_values(["assignment", "label_i", "label_j"]).reset_index(drop=True)
    return df

def read_in_tsv(file, id_to_chr):
    col_names = ['chrom',
                 'centering_position',
                 'strand',
                 'subset_sequence',
                 'reference_start',
                 'reference_end',
                 'query_name',
                 'RG',
                 'HP',
                 'centered_query_start',
                 'centered_query_end',
                 'query_length',
                 '5mC_pos',
                 '5mC_score',
                 'query_sequence']
    sample_tsv = pd.read_csv(
        file, sep="\t", index_col=False, header=None,
        usecols=[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 15, 21],
        names=col_names, skiprows=1
    )
    sample_tsv['chr_base'] = sample_tsv['chrom'].map(id_to_chr)
    threshold = 100_000
    sample_tsv['chrom_assignment'] = (
        sample_tsv['chr_base'] + np.where(sample_tsv['centering_position'].astype(int) > threshold, 'q', 'p')
    )
    sample_tsv.drop(columns=['chr_base'], inplace=True)
    return sample_tsv

if __name__ == "__main__":
    if len(sys.argv) < 3:
        raise SystemExit(
            "Usage: python script.py '<glob_pattern>' <chrom_alias.tsv> [NUM_THREADS]\n"
            "Example: python script.py 'data/*.tsv.gz' chrom_alias.tsv 16"
        )

    chrom_alias = pd.read_csv(sys.argv[2],
                              sep="\t",
                              header=None,
                              names=['id', 'chrom', 'contig'])
    id_to_chr = dict(zip(chrom_alias['id'], chrom_alias['chrom']))

    NUM_THREADS = (
        int(sys.argv[3]) if len(sys.argv) > 3
        else int(os.environ.get("THREADS", os.cpu_count() or 8))
    )
    NUM_THREADS = max(1, min(32, NUM_THREADS))
    
    glob_pattern = sys.argv[1]

    files = glob(glob_pattern)

    def _load_and_tvr(path: str) -> pd.DataFrame:
        tsv = read_in_tsv(path, id_to_chr)
        return assignTVRAndAddTelomereSequenceColumn(tsv)

    pb_hprc_tsvs = []
    with ThreadPoolExecutor(max_workers=NUM_THREADS) as ex:
        futs = {ex.submit(_load_and_tvr, f): f for f in files}
        for fut in tqdm(as_completed(futs), total=len(futs), desc=f"TSVs (threads={NUM_THREADS})"):
            pb_hprc_tsvs.append(fut.result())

    
    canonical_arms = [
        f"chr{i}{arm}" for i in range(1, 23) for arm in ("p", "q")
    ] + ["chrXp", "chrXq", "chrYp", "chrYq"]
        
        
    w_ranges = range(50, 4_050, 50)


    def process_arm(ca: str) -> pd.DataFrame:
        df_cons = collect_consensus_telomeres(
            all_tsvs=pb_hprc_tsvs,
            chrom_end=ca,
            min_reads_per_group=5,
            only_labels=None,
            t_cap=max(w_ranges)
        )
        df_sub = df_cons.copy()
        df_sub["consensus"] = df_sub["consensus"].str[:max(w_ranges)]
        pairwise = pairwise_tvr_similarity_alignment(
            df_sub, w_ranges,
            min_tvr_size=1, max_tvr_size=300, k=6
        )
        pairwise = pairwise.assign(assignment=ca).reset_index(drop=True)
        for t_cap in w_ranges:
            col = f"sim_t{t_cap}"
            pairwise[col] = pairwise["sim"].apply(lambda d: d.get(t_cap) if isinstance(d, dict) else np.nan)
        return pairwise


    pairwise_all = []
    with ProcessPoolExecutor(max_workers=NUM_THREADS) as ex:
        futs = {ex.submit(process_arm, ca): ca for ca in canonical_arms[:-2]}
        for fut in tqdm(as_completed(futs), total=len(futs), desc=f"Arms (procs={NUM_THREADS})"):
            pairwise_all.append(fut.result())

    big_pairwise = pd.concat(pairwise_all, ignore_index=True)
    keep_cols = ["assignment", "contig_i", "contig_j", "label_i", "label_j"] + [f"sim_t{t}" for t in w_ranges]
    big_pairwise = big_pairwise[keep_cols].copy()
    big_pairwise.to_csv("ont_pairwise_all_arms.csv", index=False)
    print("Wrote pairwise_all_arms.csv")

