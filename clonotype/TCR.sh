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
echo "indexing started"
samtools index -@ $threads "$sample"
sleep 5
echo "indexing done"

# get unmapped reads
echo "unmapped started"
samtools view -b -h -f 4 -@ $threads "$sample" > "$sample".unmapped.bam
sleep 5
echo "unmapped done"

# get TRA reads
echo "TRA started"
samtools view -b -h -@ $threads "$sample" "chr14:21621904-22552132" > "$sample".TRA.bam
sleep 5
echo "TRA done"

# get TRB reads
echo "TRB started"
samtools view -b -h -@ $threads "$sample" "chr7:142299011-142813287" > "$sample".TRB.bam
sleep 5
echo "TRB done"

# get TRG reads
echo "TRG started"
samtools view -b -h -@ $threads "$sample" "chr7:38240024-38368055" > "$sample".TRG.bam
sleep 5
echo "TRG done"

# get TRD reads
echo "TRG started"
samtools view -b -h -@ $threads "$sample" "chr14:22422546-22466577" > "$sample".TRD.bam
sleep 5
echo "TRG done"

# merge TR reads
echo "merge started"
samtools merge -f -@ $threads "$sample".TR.bam "$sample".TRA.bam "$sample".TRB.bam "$sample".TRG.bam "$sample".TRD.bam
sleep 5
echo "merge done"

# sort TR reads
echo "sort started"
samtools sort -@ $threads -n "$sample".TR.bam -o "$sample".TR.sorted.bam
sleep 5
echo "sort done"

echo "more processing started"

# sort unmapped reads
samtools sort -n "$sample".unmapped.bam -o "$sample".unmapped.sorted.bam

# get unmapped fastq
bedtools bamtofastq -i "$sample".unmapped.sorted.bam -fq "$sample".unmapped.R1.fastq -fq2 "$sample".unmapped.R2.fastq 2> tmp

# get TR fastq
bedtools bamtofastq -i "$sample".TR.sorted.bam -fq "$sample".TR.R1.fastq -fq2 "$sample".TR.R2.fastq 2> tmp

# get combined fastq
cat "$sample".unmapped.R1.fastq "$sample".TR.R1.fastq > "$sample".combined.T.R1.fastq
cat "$sample".unmapped.R2.fastq "$sample".TR.R2.fastq > "$sample".combined.T.R2.fastq

sleep 5
echo "more processing done"

echo "mixcr started"

echo "copy license"
cp /data/mi.license /opt/mixcr/mi.license

echo "mixcr align"
# mixcr align
mixcr align -f -r "$report" -s human -t $threads -p rna-seq -OallowPartialAlignments=true -OvParameters.geneFeatureToAlign=VGeneWithP "$sample".combined.T.R1.fastq "$sample".combined.T.R2.fastq "$sample".T.vdjca

echo "mixcr assemblePartial"
# mixcr assemblePartial
mixcr assemblePartial -f "$sample".T.vdjca "$sample".T.rescued.vdjca
mixcr assemblePartial -f "$sample".T.rescued.vdjca "$sample".T.rescued2.vdjca

echo "mixcr assemble"
# mixcr assemble
mixcr assemble --write-alignments -f "$sample".T.rescued2.vdjca "$sample".T.clna

echo "mixcr contig assembly"
mixcr assembleContigs "$sample".T.clna "$sample".T.clns

echo "mixcr exportClones"
# mixcr exportClones
mixcr exportClones -count -fraction -vGene -dGene -jGene -vAlignment -dAlignment -jAlignment -aaFeature CDR3 -p fullImputed "$sample".T.clns "$sample".T.clones.tsv
# mixcr exportClones -f -o -t "$sample".T.clns "$sample".T.clones.tsv

sleep 5
echo "mixcr done"

echo "parse output started"

# parse output
cat "$sample".T.clones.tsv | grep -v 'IGH' | grep -v 'IGK' | grep -v 'IGL' > tmp.txt

# get all outputs
head -1 tmp.txt | sed -e "s/^/sample	/" > "$all_output"
tail -n +2 tmp.txt | sed -e "s/^/$id	/" >> "$all_output"

# get header
head -1 "$all_output" > "$sample".T.clones.tsv.head1

# get top tra line
grep -m 1 "TRA" "$all_output" > "$sample".T.clones.tsv.TRA

# get top trb line 
grep -m 1 "TRB" "$all_output" > "$sample".T.clones.tsv.TRB

# get top trg line
grep -m 1 "TRG" "$all_output" > "$sample".T.clones.tsv.TRG

# get top trd line
grep -m 1 "TRD" "$all_output" > "$sample".T.clones.tsv.TRD

# cat all files together
cat "$sample".T.clones.tsv.head1 "$sample".T.clones.tsv.TRA "$sample".T.clones.tsv.TRB "$sample".T.clones.tsv.TRG "$sample".T.clones.tsv.TRD > "$top_output"

sleep 5
echo "parse output done"
