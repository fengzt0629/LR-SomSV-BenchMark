#######################
#minimap2
#! /bin/bash
fq=$1
sampleid=$2
workdir=$3
fqdir=$4
ref=$5

clean_fq=$workdir/clean.$sampleid.fq.gz
sam=$workdir/${sampleid}.sam
bam=$workdir/${sampleid}.bam
clean_fastq_plots=$workdir/${sampleid}.clean-fastq-plots
bam_sort=$workdir/${sampleid}_minimap2.2.17_sorted.bam
bam_sort_tag=$workdir/${sampleid}_minimap2.2.17_sorted_tag.bam




echo "$(date) 3. Start mapping: $sampleid"
minimap2  --MD -ax map-ont -L -t 40 $ref $fq > $sam
echo "$(date) 3. Finish mapping: $sampleid"

echo "$(date) 4. Start bam: $sampleid"
samtools view -@ 20 -bS $sam -o $bam
echo "$(date) 4. Finish bam: $sampleid"

echo "$(date) 5. Start Delete sam: $sampleid"
rm $sam
echo "$(date) 5. Finish Delete sam: $sampleid"

echo "$(date) 6. Start sort: $sampleid"
samtools sort  -@ 20 $bam -o $bam_sort_tag
echo "$(date) 6. Finish sort: $sampleid"

echo "$(date) 7. Start add tag: $sampleid"
(samtools view -H "$bam_sort"; samtools view "$bam_sort"|awk 'BEGIN{OFS="\t"}{$1="blood"$1; print $0}')|samtools view -S -b - -o "$bam_sort_tag"
echo "$(date) 7. Finish add tag: $sampleid"

echo "$(date) 8. Start Delete bam: $sampleid"
rm $bam
echo "$(date) 8. Finish Delete bam: $sampleid"

echo "$(date) 9.start index sort tag bam: $sampleid"
samtools index  $bam_sort_tag
echo "$(date) 9. end index sort tag bam : $sampleid"
