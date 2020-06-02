old_bam=$1
threads=$2
flag=$3
new_bam=$4
samtools view -h -F ${flag} ${old_bam} | awk '$7!="="' | samtools view -b -@ ${threads} -S - > ${new_bam}
samtools index ${new_bam}