# the point of this script is to get mpileup filter flags
# from cleaned variant filtering outputs

bam=$1
clean_vars=$2
mpileup_vars=$3
sample_id=$4
depth_thresh=$5
std_thresh=$6
ref=$7

# remove headers and get column names
grep -v 'CHROM_POS_REF_ALT' $clean_vars | cut -f2,3 > all_calls.txt

# run samtools mpileup
samtools mpileup --output-BP -a -f $ref -l all_calls.txt $bam > all_vars.txt

python3 add_mpileup_af_20211202.py all_vars.txt $clean_vars $depth_thresh $std_thresh $mpileup_vars