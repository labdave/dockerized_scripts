while getopts ":w:d:g:s:o:" opt; do
  case $opt in
    w) whitelist="$OPTARG"
    ;;
    d) dv="$OPTARG"
    ;;
    g) hc="$OPTARG"
    ;;
    s) st="$OPTARG"
    ;;
    o) output="$OPTARG"
	;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

### DEEPVARIANT ###
gunzip -f "${dv}".gz
sed '/##/d' "${dv}" > "${dv}"_dv_noheader.vcf
python reformat_deepvariant.py "${dv}"_dv_noheader.vcf "${dv}"_dv_reformat.vcf
python intersect_vcf.py "${whitelist}" "${dv}"_dv_reformat.vcf "${dv}"_dv_intersect.vcf
echo "DV done"

### STRELKA2 ###
gunzip -f "${st}".gz
sed '/##/d' "${st}" > "${st}"_s2_noheader.vcf
python reformat_strelka.py "${st}"_s2_noheader.vcf "${st}"_s2_reformat.vcf
python intersect_vcf.py "${whitelist}" "${st}"_s2_reformat.vcf "${st}"_s2_intersect.vcf
echo "S2 done"

### HAPLOTYPECALLER ###
gunzip -f "${hc}".gz
sed '/##/d' "${hc}" > "${hc}"_hc_noheader.vcf
python reformat_haplotypecaller.py "${hc}"_hc_noheader.vcf "${hc}"_hc_reformat.vcf
python intersect_vcf.py "${whitelist}" "${hc}"_hc_reformat.vcf "${hc}"_hc_intersect.vcf
echo "HC done"

### MERGER ###
python final_parser.py "${dv}"_dv_intersect.vcf "${st}"_s2_intersect.vcf "${hc}"_hc_intersect.vcf "${output}"