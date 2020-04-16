# obtain parameters
bam=$1
sample=$(echo "${bam}" | sed 's/.bam//g')
read_length=$2
discordant_z=$3
back_distance=$4
weight=$5
min_mapping_threshold=$6

# locate lumpy scripts directory
scripts="lumpy-sv/scripts"

# index bam file
samtools index "${sample}".bam

# extract the discordant paired-end alignments.
samtools view -b -F 1294 "${sample}".bam > "${sample}".discordant.bam

# extract the split-read alignments
samtools view -h "${sample}".bam | "${scripts}"/extractSplitReads_BwaMem -i stdin | samtools view -Sb - > "${sample}".split.bam

# sort both alignments
samtools sort "${sample}".discordant.bam > "${sample}".discordant.sorted.bam
samtools index "${sample}".discordant.sorted.bam

samtools sort "${sample}".split.bam > "${sample}".split.sorted.bam
samtools index "${sample}".split.sorted.bam

# generate empirical insert size statistics on each library in the BAM file 
samtools view "${sample}".bam | tail -n+100000 | "${scripts}"/pairend_distro.py -r "${read_length}" -X 4 -N 10000 -o "${sample}".lib.histo > tmp.out 2>&1
mean=`tail -1 tmp.out | cut -f1 | cut -f2 -d ':' | cut -f1 -d '.'`
std=`tail -1 tmp.out | cut -f2 | cut -f2 -d ':' | cut -f1 -d '.'`

# lumpy call (traditional)
lumpy -mw 4 -tt 0 \
-pe id:sample,bam_file:"${sample}".discordant.sorted.bam,histo_file:"${sample}".lib.histo,mean:$mean,stdev:$std,read_length:"${read_length}",min_non_overlap:"${read_length}",discordant_z:"${discordant_z}",back_distance:"${back_distance}",weight:"${weight}",min_mapping_threshold:"${min_mapping_threshold}" \
-sr id:sample,bam_file:"${sample}".split.sorted.bam,back_distance:"${back_distance}",weight:"${weight}",min_mapping_threshold:"${min_mapping_threshold}" \
> "${sample}".lumpy.vcf

# post procesing with SVTyper to make GT Calls from the Lumpy vcf  using a Bayesian maximum likelihood algorithm.
svtyper -B "${sample}".bam -S "${sample}".split.sorted.bam -i "${sample}".lumpy.vcf > "${sample}".gt.vcf