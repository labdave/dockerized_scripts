# remove header
awk -v FS='	' -v OFS='	' '{print $2,$3,$4,$5,$6,$1}' $1 | sed -i '/chrom/d' > $1.swp

# get complement
bedtools complement -i $1.swp -g hg38.fa.fai > $1.comp

# add 0 cnv and sample id for complement regions
awk -v FS='	' -v OFS='	' -v ref=$3 '{print $1,$2+1,$3-1,0,0,ref}' $1.comp > tmp

# combine both files
cat tmp $1.swp > tmp1

# sort combined file
bedtools sort -i tmp1 > tmp

# move id to first column
awk -v FS='	' -v OFS='	' '{print $6,$1,$2,$3,$4,$5}' tmp > $2