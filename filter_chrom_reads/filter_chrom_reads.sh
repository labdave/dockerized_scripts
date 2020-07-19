old_bam=$1
threads=$2
flag=$3
new_bam=$4
temp_bam1=$5
temp_bam2=$6
rm ${new_bam}
samtools view -@ ${threads} -h -F ${flag} ${old_bam} | awk '$7!="="' | samtools view -b -@ ${threads} -S - > ${temp_bam1}
samtools view -@ ${threads} -h ${old_bam} | head -100000 | samtools view -b -@ ${threads} -S - > ${temp_bam2}

# add non-primary alignments
samtools view -@ ${threads} -h -f 256 ${old_bam} | samtools view -Sb -@ ${threads} - > non-primary.bam 

samtools merge -@ ${threads} ${new_bam} ${temp_bam1} ${temp_bam2} non-primary.bam
samtools index -@ ${threads} ${new_bam}