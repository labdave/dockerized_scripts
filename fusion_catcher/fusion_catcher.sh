input_dir=$1
output_dir=$2
data_dir=$3

/opt/fusioncatcher/bin/fusioncatcher.py -i "${input_dir}" -o "${output_dir}" -d "${data_dir}"
