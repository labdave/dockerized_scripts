#!/bin/bash
old_bam=$1
threads=$2
flag=$3
new_bam=$4
padded_bed=$5
rm ${new_bam}

echo "start getting header"
time samtools view -H ${old_bam} > /data/output/header.sam
echo "end getting header"

echo "start filtering by flag (-F, exclude ${flag})"
time samtools view -@ ${threads} -F ${flag} ${old_bam} | cut -d'	' -f1 > /data/output/tmp.flag_reads
echo "end filtering by flag"

echo "start filtering for good reads (-f 2, include first 100,000 properly paired reads)"
time samtools view -@ ${threads} -f 2 ${old_bam} | cut -d'	' -f1 | head -100000 > /data/output/tmp.good_reads
echo "end filtering for good reads"

echo "start filtering for non-primary reads (-f 256)"
time samtools view -@ ${threads} -f 256 ${old_bam} | cut -d'	' -f1 > /data/output/tmp.non_primary
echo "end filtering for non-primary reads"

echo "combine read IDs from all filtering steps"
cat /data/output/tmp.flag_reads /data/output/tmp.good_reads /data/output/tmp.non_primary > /data/output/tmp.all_read_ids.txt
# sort and get unique read IDs
time sort -u --parallel ${threads} /data/output/tmp.all_read_ids.txt > /data/output/tmp.unique_read_ids.txt
# split into multiple files, at most 5,000,000 records per file
time split -l 5000000 --numeric-suffixes /data/output/tmp.unique_read_ids.txt /data/output/split_2_files

echo "start subsetting BAM file by read ID file(s)"
FILES=/data/output/split_2_files*
for filename in ${FILES}; do
	echo $filename
	name=$( echo ${filename} | cut -d'/' -f4 )
	echo $name
	time java -jar picard.jar FilterSamReads I=${old_bam} O=/data/output/bam2_${name}.bam RLF=${filename} FILTER=includeReadList &
done
wait $(jobs -p)
echo "end subsetting BAM file by read ID file(s)"

echo "start merging"
time samtools merge -fcp -@ ${threads} -h /data/output/header.sam ${new_bam} /data/output/bam2_*
echo "end merging"

echo "start removing hard clipped reads"
samtools view -@ ${threads} ${new_bam} | awk '$6 !~ /.*H.*H/ {print $0}' > /data/output/tmp.sam
cat /data/output/header.sam /data/output/tmp.sam | samtools view -Sb - > ${new_bam}

echo "start indexing"
time samtools index -@ ${threads} ${new_bam}
echo "end indexing"
