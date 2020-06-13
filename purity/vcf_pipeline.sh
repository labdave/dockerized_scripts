ext_flag="FALSE"
while getopts ":w:e:d:g:s:i:" opt; do
  case $opt in
    w) whitelist="$OPTARG"
    ;;
    e) extended_whitelist="$OPTARG"; ext_flag="TRUE"
    ;;
    d) dv="$OPTARG"
    ;;
    g) hc="$OPTARG"
    ;;
    s) st="$OPTARG"
    ;;
    i) sample_id="$OPTARG"
	;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

### DEEPVARIANT ###
gunzip -f "${dv}".gz
sed '/##/d' "${dv}" > "${dv}"_dv_noheader.vcf
python reformat_deepvariant.py "${dv}"_dv_noheader.vcf "${dv}"_dv_reformat.vcf
if ((ext_flag == "TRUE"))
then
	python intersect_vcf.py "${whitelist}" "${dv}"_dv_reformat.vcf "${dv}"_dv_intersect.vcf "${extended_whitelist}"
else
	python intersect_vcf.py "${whitelist}" "${dv}"_dv_reformat.vcf "${dv}"_dv_intersect.vcf
fi

echo "DV done"

### STRELKA2 ###
gunzip -f "${st}".gz
sed '/##/d' "${st}" > "${st}"_s2_noheader.vcf
python reformat_strelka.py "${st}"_s2_noheader.vcf "${st}"_s2_reformat.vcf
if ((ext_flag == "TRUE"))
then
	python intersect_vcf.py "${whitelist}" "${st}"_s2_reformat.vcf "${st}"_s2_intersect.vcf "${extended_whitelist}"
else
	python intersect_vcf.py "${whitelist}" "${st}"_s2_reformat.vcf "${st}"_s2_intersect.vcf
fi

echo "S2 done"

### HAPLOTYPECALLER ###
gunzip -f "${hc}".gz
sed '/##/d' "${hc}" > "${hc}"_hc_noheader.vcf
python reformat_haplotypecaller.py "${hc}"_hc_noheader.vcf "${hc}"_hc_reformat.vcf
if ((ext_flag == "TRUE"))
then
	python intersect_vcf.py "${whitelist}" "${hc}"_hc_reformat.vcf "${hc}"_hc_intersect.vcf "${extended_whitelist}"
else
	python intersect_vcf.py "${whitelist}" "${hc}"_hc_reformat.vcf "${hc}"_hc_intersect.vcf
fi

echo "HC done"

### MERGER ###
python final_parser.py "${dv}"_hc_intersect.vcf "${s2}"_hc_intersect.vcf "${hc}"_hc_intersect.vcf "${sample_id}".purity.tsv