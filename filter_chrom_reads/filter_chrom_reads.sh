old_bam=$1
threads=$2
flag=$3
new_bam=$4
temp_bam1=$5
temp_bam2=$6
rm ${new_bam}
samtools view -@ ${threads} -h -F ${flag} ${old_bam} | awk '$7!="="' | samtools view -b -@ ${threads} -S - > ${temp_bam1}
samtools view -@ ${threads} -h ${old_bam} | head -10000 | samtools view -b -@ ${threads} -S - > ${temp_bam2}
samtools merge -@ ${threads} ${new_bam} ${temp_bam1} ${temp_bam2}
samtools index -@ ${threads} ${new_bam}