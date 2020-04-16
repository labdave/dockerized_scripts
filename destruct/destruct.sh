# obtain parameters
bam=$1
cfg=$2

destruct_ref_data=/home/naresh_prodduturi/data/HG38/destruct_ref/
rundir=/home/naresh_prodduturi/data/translocations/Date_1_3_2020_8MB/Delly_ALT/destruct
echo "$rundir/$i"
mkdir -p $rundir/$i
cd $rundir/$i
input_bam=/home/naresh_prodduturi/data/translocations/Date_1_3_2020_8MB/Delly_ALT/bamfiles/$i.bwamem_alt_snaut.bam
destruct run  $destruct_ref_data breaks.tsv break_libs.tsv break_reads.tsv --bam_files $input_bam --lib_ids sample.$i --raw_data_dir destruct_raw_data/ --submit local --config $cfg