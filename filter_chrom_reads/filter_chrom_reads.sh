old_bam=$1
threads=$2
flag=$3
new_bam=$4
temp_bam1=$5
temp_bam2=$6
padded_bed=$7
rm ${new_bam}

echo "started filtering for on-target reads"
time samtools view -@ ${threads} -L ${padded_bed} ${old_bam} | cut -d'	' -f1 > tmp
time sort -u -S20G --parallel ${threads} tmp > on-target.reads.txt
time java -jar picard.jar FilterSamReads I=${old_bam} O=on-target.bam RLF=on-target.reads.txt FILTER=includeReadList
echo "ended filtering for on-target reads"

echo "started filtering by flag"
time samtools view -@ ${threads} -h -F ${flag} ${on-target.bam} | awk '$7!="="' | samtools view -b -@ ${threads} -S - > ${temp_bam1}
echo "ended filtering by flag"

echo "started filtering for good reads"
time samtools view -@ ${threads} -h on-target.bam | head -100000 | samtools view -b -@ ${threads} -S - > ${temp_bam2}
echo "ended filtering for good reads"

echo "started filtering for non-primary reads"
# add non-primary alignments
time samtools view -@ ${threads} -f 256 on-target.bam | cut -d'	' -f1 > tmp
time sort -u -S20G --parallel ${threads} tmp > non-primary.reads.txt
time java -jar picard.jar FilterSamReads I=on-target.bam O=non-primary.bam RLF=non-primary.reads.txt FILTER=includeReadList
echo "ended filtering for non-primary reads"

echo "started merging"
time samtools merge -cp -@ ${threads} ${new_bam} ${temp_bam1} ${temp_bam2} non-primary.bam
echo "ended merging"

echo "started indexing"
time samtools index -@ ${threads} ${new_bam}
echo "ended indexing"