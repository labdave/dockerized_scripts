#!/bin/bash
old_bam=$1
threads=$2
flag=$3
new_bam=$4
padded_bed=$5
rm ${new_bam}

echo "started getting header"
time samtools view -H ${old_bam} > /data/output/header.sam
echo "ended getting header"

echo "started filtering by flag"
time samtools view -@ ${threads} -h -F ${flag} ${old_bam} | awk '$7!="="' | samtools view -b -@ ${threads} -S - > /data/output/temp1.bam
echo "ended filtering by flag"

echo "started filtering for good reads"
time samtools view -@ ${threads} -h ${old_bam} | head -100000 | samtools view -b -@ ${threads} -S - > /data/output/temp2.bam
echo "ended filtering for good reads"

echo "started filtering for non-primary reads"
# add non-primary alignments
time samtools view -@ ${threads} -f 256 ${old_bam} | cut -d'	' -f1 > /data/output/tmp1
time sort -u -S20G --parallel ${threads} /data/output/tmp1 > /data/output/non-primary.reads.txt
time split -l 5000000 --numeric-suffixes /data/output/non-primary.reads.txt /data/output/split_2_files
FILES=/data/output/split_2_files*
for filename in ${FILES}; do
	echo $filename
	name=$( echo ${filename} | cut -d'/' -f4 )
	echo $name
	time java -jar picard.jar FilterSamReads I=${old_bam} O=/data/output/bam2_${name}.bam RLF=${filename} FILTER=includeReadList &
done
wait $(jobs -p)
time samtools merge -fcp -@ ${threads} -h /data/output/header.sam /data/output/non-primary.bam /data/output/bam2_*
echo "ended filtering for non-primary reads"

echo "started merging"
time samtools merge -fcp -@ ${threads} ${new_bam} /data/output/temp1.bam /data/output/temp2.bam /data/output/non-primary.bam
echo "ended merging"

echo "start removing hard clipped reads"
samtools view -@ ${threads} ${new_bam} | awk 'match($6, !/.*H.*H/) {print $0}' | samtools view -h -Sb - /data/output/header.sam > /data/output/tmp
samtools merge -fcp -@ ${threads} -h /data/output/header.sam ${new_bam} /data/output/tmp

echo "started indexing"
time samtools index -@ ${threads} ${new_bam}
echo "ended indexing"