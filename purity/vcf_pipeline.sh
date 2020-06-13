ext_flag="FALSE"
while getopts ":w:e:d:s:h" opt; do
  case $opt in
    w) whitelist="$OPTARG"
    ;;
    e) extended_whitelist="$OPTARG"; ext_flag="TRUE"
    ;;
    d) dv="$OPTARG"
    ;;
    s) s2="$OPTARG"
    ;;
    h) hc="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

### DEEPVARIANT ###
gunzip -f "${dv}"
sed '/##/d' "${dv}" > "${dv}"_dv_noheader.vcf
python reformat_deepvariant.py "${dv}"_dv_noheader.vcf "${dv}"_dv_reformat.vcf
if ((ext_flag == "TRUE"))
then
	python intersect_vcf.py "${whitelist}" "${dv}"_dv_reformat.vcf "${sample_id}"_dv_intersect.vcf "${extended_whitelist}"
else
	python intersect_vcf.py "${whitelist}" "${dv}"_dv_reformat.vcf "${sample_id}"_dv_intersect.vcf
fi

### STRELKA2 ###
gunzip -f "${s2}"
sed '/##/d' "${s2}" > "${s2}"_s2_noheader.vcf
python reformat_strelka.py "${s2}"_s2_noheader.vcf "${s2}"_s2_reformat.vcf
if ((ext_flag == "TRUE"))
then
	python intersect_vcf.py "${whitelist}" "${s2}"_s2_reformat.vcf "${sample_id}"_s2_intersect.vcf "${extended_whitelist}"
else
	python intersect_vcf.py "${whitelist}" "${s2}"_s2_reformat.vcf "${sample_id}"_s2_intersect.vcf
fi

### HAPLOTYPECALLER ###
gunzip -f "${hc}"
sed '/##/d' "${hc}" > "${hc}"_hc_noheader.vcf
python reformat_haplotypecaller.py "${hc}"_hc_noheader.vcf "${hc}"_hc_reformat.vcf
if ((ext_flag == "TRUE"))
then
	python intersect_vcf.py "${whitelist}" "${hc}"_hc_reformat.vcf "${sample_id}"_hc_intersect.vcf "${extended_whitelist}"
else
	python intersect_vcf.py "${whitelist}" "${hc}"_hc_reformat.vcf "${sample_id}"_hc_intersect.vcf
fi