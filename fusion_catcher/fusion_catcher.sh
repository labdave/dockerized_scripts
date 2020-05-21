input_dir=$1
output_dir=$2
data_dir="/opt/fusioncatcher/data/human_v98/"

/opt/fusioncatcher/bin/fusioncatcher.py -i "${input_dir}" -o "${output_dir}" -d "${data_dir}"
