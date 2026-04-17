#!/usr/bin/env bash
set -euo pipefail

# Compare NM (edit distance) and AS (alignment score) for:
#   1) telomere-spanning reads  (reads overlapping telomere BED)
#   2) non-telomere-spanning reads (all other mapped reads)
#
# Usage:
#   ./telo_nm_as.sh reads.bam telomeres.bed > nm_as_summary.tsv
#
# Notes:
# - Requires: samtools, awk
# - BAM should be indexed (script will try to index if missing)
# - BED should match BAM contig names

BAM="${1:?need BAM}"
BED="${2:?need telomere BED}"
THREADS="${THREADS:-8}"

if [[ ! -f "${BAM}.bai" && ! -f "${BAM%.bam}.bai" ]]; then
  samtools index -@ "${THREADS}" "${BAM}"
fi

tmpdir="$(mktemp -d)"
trap 'rm -rf "${tmpdir}"' EXIT

telo_names="${tmpdir}/telo.qnames.txt"
all_names="${tmpdir}/all.qnames.txt"
nontelo_names="${tmpdir}/nontelo.qnames.txt"

# Telomere-overlapping mapped read names (unique)
samtools view -@ "${THREADS}" -F 4 -L "${BED}" "${BAM}" \
  | awk '{print $1}' | sort -u > "${telo_names}"

# All mapped read names (unique)
samtools view -@ "${THREADS}" -F 4 "${BAM}" \
  | awk '{print $1}' | sort -u > "${all_names}"

# Non-telomere mapped read names (set difference)
comm -23 "${all_names}" "${telo_names}" > "${nontelo_names}"

# Helper: stream alignments for qname list and summarize NM/AS
summarize_nm_as () {
  local label="$1"
  local qnames="$2"

  # samtools view can filter by read-name list with -N
  samtools view -@ "${THREADS}" -F 4 -N "${qnames}" "${BAM}" \
  | awk -v LABEL="${label}" '
    function get_tag(prefix,   i, v) {
      # tags begin at field 12 in SAM
      for (i=12; i<=NF; i++) {
        if (index($i, prefix)==1) {
          v = substr($i, length(prefix)+1)
          return v
        }
      }
      return ""
    }
    BEGIN {
      n=0;
      nm_n=0; nm_sum=0; nm_sum2=0;
      as_n=0; as_sum=0; as_sum2=0;
    }
    {
      n++;

      nm = get_tag("NM:i:");
      if (nm != "") {
        nm_n++; nm_sum += nm; nm_sum2 += nm*nm;
      }

      as = get_tag("AS:i:");
      if (as != "") {
        as_n++; as_sum += as; as_sum2 += as*as;
      }
    }
    END {
      nm_mean = (nm_n>0 ? nm_sum/nm_n : "NA");
      nm_sd   = (nm_n>1 ? sqrt((nm_sum2/nm_n) - (nm_sum/nm_n)^2) : "NA");
      as_mean = (as_n>0 ? as_sum/as_n : "NA");
      as_sd   = (as_n>1 ? sqrt((as_sum2/as_n) - (as_sum/as_n)^2) : "NA");

      # Output a single TSV row
      # label, N_alignments, NM_n, NM_mean, NM_sd, AS_n, AS_mean, AS_sd
      printf("%s\t%d\t%d\t%s\t%s\t%d\t%s\t%s\n",
             LABEL, n, nm_n, nm_mean, nm_sd, as_n, as_mean, as_sd);
    }'
}

# Header
echo -e "group\tN_alignments\tNM_n\tNM_mean\tNM_sd\tAS_n\tAS_mean\tAS_sd"

summarize_nm_as "telomere_overlapping" "${telo_names}"
summarize_nm_as "non_telomere_overlapping" "${nontelo_names}"

