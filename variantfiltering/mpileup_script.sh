# the point of this script is to get mpileup filter flags
# from cleaned variant filtering outputs

bam=$1
clean_vars=$2
oncokb_genes=$3
mpileup_vars=$4
sample_id=$5
depth_thresh=$6
std_thresh=$7
ref=$8

# remove headers and get column names
grep -v 'CHROM_POS_REF_ALT' $clean_vars | cut -f2,3 > all_calls.txt

# run samtools mpileup
samtools mpileup --output-BP -a -f $ref -l all_calls.txt $bam > all_vars.txt

python3 add_mpileup_af.py all_vars.txt $clean_vars $depth_thresh $std_thresh $mpileup_vars