# get the areas of interest from bam file into multiple bams based on conditions
# get mates somehow and make into paired end fastq

sample=$1
threads=$2

report=$3
top_output=$4
all_output=$5
id=$6

echo "ARGUMENTS"
echo "sample:" "$sample"
echo "threads:" "$threads"
echo "report:" "$report"
echo "top_output:" "$top_output"
echo "all_output:" "$all_output"
echo "sample id:" "$id"
echo "========="
echo

# get index
#samtools index -@ $threads "$sample"

# get unmapped reads
samtools view -b -h -f 4 -@ $threads "$sample" > "$sample".unmapped.bam

# get TRA reads
samtools view -b -h -@ $threads "$sample" "chr14:21621904-22552132" > "$sample".TRA.bam

# get TRB reads
samtools view -b -h -@ $threads "$sample" "chr7:142299011-142813287" > "$sample".TRB.bam

# get TRG reads
samtools view -b -h -@ $threads "$sample" "chr7:38240024-38368055" > "$sample".TRG.bam

# get TRD reads
samtools view -b -h -@ $threads "$sample" "chr14:22422546-22466577" > "$sample".TRD.bam

# merge TR reads
samtools merge -f -@ $threads "$sample".TR.bam "$sample".TRA.bam "$sample".TRB.bam "$sample".TRG.bam "$sample".TRD.bam

# sort TR reads
samtools sort -@ $threads -n "$sample".TR.bam "$sample".TR.sorted

# sort unmapped reads
samtools sort -n "$sample".unmapped.bam "$sample".unmapped.sorted

# get unmapped fastq
bedtools bamtofastq -i "$sample".unmapped.sorted.bam -fq "$sample".unmapped.R1.fastq -fq2 "$sample".unmapped.R2.fastq 2> tmp

# get TR fastq
bedtools bamtofastq -i "$sample".TR.sorted.bam -fq "$sample".TR.R1.fastq -fq2 "$sample".TR.R2.fastq 2> tmp

# get combined fastq
cat "$sample".unmapped.R1.fastq "$sample".TR.R1.fastq > "$sample".combined.T.R1.fastq
cat "$sample".unmapped.R2.fastq "$sample".TR.R2.fastq > "$sample".combined.T.R2.fastq

# mixcr align
mixcr align -f -r "$report" -s human -t $threads -p rna-seq -OallowPartialAlignments=true -OvParameters.geneFeatureToAlign=VGeneWithP "$sample".combined.T.R1.fastq "$sample".combined.T.R2.fastq "$sample".T.vdjca

# mixcr assemblePartial
mixcr assemblePartial -f "$sample".T.vdjca "$sample".T.rescued.vdjca
mixcr assemblePartial -f "$sample".T.rescued.vdjca "$sample".T.rescued2.vdjca

# mixcr assemble
mixcr assemble -f "$sample".T.rescued2.vdjca "$sample".T.clns

# mixcr exportClones
mixcr exportClones -count -fraction -vGene -dGene -jGene -vAlignment -dAlignment -jAlignment -aaFeature CDR3 "$sample".T.clns "$sample".T.clones.tsv
# mixcr exportClones -f -o -t "$sample".T.clns "$sample".T.clones.tsv

# parse output
cat "$sample".T.clones.tsv | grep -v 'IGH' | grep -v 'IGK' | grep -v 'IGL' > "$all_output"
# get header
head -1 "$all_output" > "$sample".T.clones.tsv.head1
# add sample info header
sed -i -e "s/^/sample	/" "$sample".T.clones.tsv.head1
# get top tra line
grep -m 1 "TRA" "$all_output" > "$sample".T.clones.tsv.TRA
# add sample info
sed -i -e "s/^/$id	/" "$sample".T.clones.tsv.TRA
# get top trb line 
grep -m 1 "TRB" "$all_output" > "$sample".T.clones.tsv.TRB
# add sample info
sed -i -e "s/^/$id	/" "$sample".T.clones.tsv.TRB
# get top trg line
grep -m 1 "TRG" "$all_output" > "$sample".T.clones.tsv.TRG
# add sample info
sed -i -e "s/^/$id	/" "$sample".T.clones.tsv.TRG
# get top trd line
grep -m 1 "TRD" "$all_output" > "$sample".T.clones.tsv.TRD
# add sample info
sed -i -e "s/^/$id	/" "$sample".T.clones.tsv.TRD
# cat all files together
grep -m 1 "TRA" "$all_output" > "$sample".T.clones.tsv.TRA
grep -m 1 "TRB" "$all_output" > "$sample".T.clones.tsv.TRB
grep -m 1 "TRG" "$all_output" > "$sample".T.clones.tsv.TRG
grep -m 1 "TRD" "$all_output" > "$sample".T.clones.tsv.TRD
cat "$sample".T.clones.tsv.head1 "$sample".T.clones.tsv.TRA "$sample".T.clones.tsv.TRB "$sample".T.clones.tsv.TRG "$sample".T.clones.tsv.TRD > "$top_output"
# cut -f 2,3,4,6,7,8,9 "$sample".T.clones.tsv.filter2 > "$sample".T.clones.tsv.filter3
# sed 's/[(][^)]*[)]//g' "$sample".T.clones.tsv.filter3 > "$output"