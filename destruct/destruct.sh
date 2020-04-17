# obtain parameters
bam=$1
sample=$(echo "${bam}" | sed 's/.bam//g')
cfg=$2
out_breaks=$3
out_break_libs=$4
out_break_reads=$5

destruct_ref_data="/usr/local/bin/destruct_ref/"
destruct run "${destruct_ref_data}" "${out_breaks}" "${out_break_libs}" "${out_break_reads}" --bam_files "${bam}" --lib_ids "${sample}" --raw_data_dir "${destruct_raw_data}" --submit local