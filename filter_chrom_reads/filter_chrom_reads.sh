old_bam=$1
threads=$2
flag=$3
new_bam=$4
padded_bed=$5
rm ${new_bam}

echo "started filtering for on_target reads"
time samtools view -H ${old_bam} > /data/output/header.sam
time samtools view -@ ${threads} -L ${padded_bed} ${old_bam} | cut -d'	' -f1 > /data/output/tmp
time sort -u -S20G --parallel ${threads} /data/output/tmp > /data/output/on_target.reads.txt
time split -l 5000000 --numeric-suffixes /data/output/on_target.reads.txt split_1_files
for filename in "/data/output/split_1_files*"; do
	name=$( echo ${filename} | cut -d'/' -f3 )
	time java -jar picard.jar FilterSamReads I=${old_bam} O=/data/output/bam_${name}.bam RLF=${filename} FILTER=includeReadList;
done
time samtools merge -cp -h /data/output/header.sam /data/output/on_target.bam /data/output/bam_*
echo "ended filtering for on_target reads"

echo "started filtering by flag"
time samtools view -@ ${threads} -h -F ${flag} /data/output/on_target.bam | awk '$7!="="' | samtools view -b -@ ${threads} -S - > /data/output/temp1.bam
echo "ended filtering by flag"

echo "started filtering for good reads"
time samtools view -@ ${threads} -h /data/output/on_target.bam | head -100000 | samtools view -b -@ ${threads} -S - > /data/output/temp2.bam
echo "ended filtering for good reads"

echo "started filtering for non-primary reads"
# add non-primary alignments
time samtools view -@ ${threads} -f 256 /data/output/on_target.bam | cut -d'	' -f1 > /data/output/tmp1
time sort -u -S20G --parallel ${threads} /data/output/tmp1 > /data/output/non-primary.reads.txt
time java -jar picard.jar FilterSamReads I=/data/output/on_target.bam O=/data/output/non-primary.bam RLF=/data/output/non-primary.reads.txt FILTER=includeReadList
echo "ended filtering for non-primary reads"

echo "started merging"
time samtools merge -cp -@ ${threads} ${new_bam} /data/output/temp1.bam /data/output/temp2.bam /data/output/non-primary.bam
echo "ended merging"

echo "started indexing"
time samtools index -@ ${threads} ${new_bam}
echo "ended indexing"