old_bam=$1
threads=$2
flag=$3
new_bam=$4
temp_bam1=$5
temp_bam2=$6
samtools view -h -F ${flag} ${old_bam} | awk '$7!="="' | samtools view -b -@ ${threads} -S - > ${temp_bam1}
samtools view -h ${old_bam} | head -1000 | samtools view -b -@ ${threads} -S - > ${temp_bam2}
samtools merge ${new_bam} ${temp_bam1} ${temp_bam2}
samtools index ${new_bam}