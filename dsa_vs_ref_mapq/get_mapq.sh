bam=$1
q=$2

samtools view -F 2308 $bam  | awk -v min_q=$q '{total++; if($5 >= min_q) pass++} END {print "Total Primary Mapped: " total+0 " | MAPQ > " min_q ": " pass+0}'
