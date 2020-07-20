old_bam=$1
threads=$2
flag=$3
new_bam=$4
temp_bam1=$5
temp_bam2=$6
rm ${new_bam}
echo "started filtering by flag"
samtools view -@ ${threads} -h -F ${flag} ${old_bam} | awk '$7!="="' | samtools view -b -@ ${threads} -S - > ${temp_bam1}
echo "ended filtering by flag"

echo "started filtering for good reads"
samtools view -@ ${threads} -h ${old_bam} | head -100000 | samtools view -b -@ ${threads} -S - > ${temp_bam2}
echo "ended filtering for good reads"

echo "started filtering for non-primary reads"
# add non-primary alignments
samtools view -@ ${threads} -f 256 ${old_bam} | cut -d'     ' -f1 | sort | uniq > non-primary.reads.txt
java -jar /home/picard.jar FilterSamReads I=${old_bam} O=non-primary.bam RLF=non-primary.reads.txt FILTER=includeReadList
echo "ended filtering for non-primary reads"

echo "started merging"
samtools merge -cp -@ ${threads} ${new_bam} ${temp_bam1} ${temp_bam2} non-primary.bam
echo "ended merging"

echo "started indexing"
samtools index -@ ${threads} ${new_bam}
echo "ended indexing"