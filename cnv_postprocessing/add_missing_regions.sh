# get complement
bedtools complement -i $1 -g hg38.fa.fai > $1.comp

# add 0 cnv and sample id for complement regions
awk -v FS='	' -v OFS='	' '{print $1,$2+1,$3-1,0,0,FILENAME}' $1.comp | sed 's/.txt.comp//g' > tmp

# combine both files
cat tmp $1 > tmp1

# sort combined file
bedtools sort -i tmp1 > tmp

# move id to first column
awk -v FS='	' -v OFS='	' '{print $6,$1,$2,$3,$4,$5}' tmp > $2