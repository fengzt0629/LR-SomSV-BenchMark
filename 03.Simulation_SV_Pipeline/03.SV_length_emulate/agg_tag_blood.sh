#minimap2
#! /bin/bash

workdir=$1

bam=$workdir/sim.srt.bam
bam_sort=$workdir/sim.srt.sort.bam
bam_sort_tag=$workdir/sim.srt.sort.tag.bam

echo "$(date) 6. Start sort: $sampleid"
samtools sort  -@ 20 $bam -o $bam_sort
echo "$(date) 6. Finish sort: $sampleid"


echo "$(date) 7. Start add tag: $sampleid"
(samtools view -H "$bam_sort"; samtools view "$bam_sort"|awk 'BEGIN{OFS="\t"}{$1="blood"$1; print $0}')|samtools view -S -b - -o "$bam_sort_tag"
echo "$(date) 7. Finish add tag: $sampleid"


echo "$(date) 6. Start index : $sampleid"
samtools index $bam_sort_tag
echo "$(date) 6. Finish idex: $sampleid"
