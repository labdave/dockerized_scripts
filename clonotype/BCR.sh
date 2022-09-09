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

# get IGH reads
echo "IGH started"
samtools view -b -h -@ $threads "$sample" "chr14:105586437-106879844" > "$sample".IGH.bam
sleep 5
echo "IGH done"

# get IGK reads
echo " IGK started"
samtools view -b -h -@ $threads "$sample" "chr2:88857361-90235368" > "$sample".IGK.bam
sleep 5
echo "IGK done"

# get IGL reads
echo "IGL started"
samtools view -b -h -@ $threads "$sample" "chr22:22026076-22922913" > "$sample".IGL.bam
sleep 5
echo "IGL done"

# merge IG reads
echo "merge started"
samtools merge -f -@ $threads "$sample".IG.bam "$sample".IGH.bam "$sample".IGK.bam "$sample".IGL.bam
sleep 5
echo "merge done"

# sort IG reads
echo "sort started"
samtools sort -@ $threads -n "$sample".IG.bam -o "$sample".IG.sorted
sleep 5
echo "sort done"

echo "more processing started"

# sort unmapped reads
samtools sort -@ $threads -n "$sample".unmapped.bam -o "$sample".unmapped.sorted

# get unmapped fastq
bedtools bamtofastq -i "$sample".unmapped.sorted.bam -fq "$sample".unmapped.R1.fastq -fq2 "$sample".unmapped.R2.fastq 2> tmp

# get IG fastq
bedtools bamtofastq -i "$sample".IG.sorted.bam -fq "$sample".IG.R1.fastq -fq2 "$sample".IG.R2.fastq 2> tmp

# debug
ls -l

# get combined fastq
cat "$sample".unmapped.R1.fastq "$sample".IG.R1.fastq > "$sample".combined.IG.R1.fastq
cat "$sample".unmapped.R2.fastq "$sample".IG.R2.fastq > "$sample".combined.IG.R2.fastq

sleep 5
echo "more processing done"

echo "mixcr started"

echo "mixcr align"
# mixcr align
mixcr align -f -r "$report" -s human -t $threads -p rna-seq -OallowPartialAlignments=true -OvParameters.geneFeatureToAlign=VGeneWithP "$sample".combined.IG.R1.fastq "$sample".combined.IG.R2.fastq "$sample".IG.vdjca

echo "mixcr assemblePartial"
# mixcr assemblePartial
mixcr assemblePartial -f "$sample".IG.vdjca "$sample".IG.rescued.vdjca
mixcr assemblePartial -f "$sample".IG.rescued.vdjca "$sample".IG.rescued2.vdjca

echo "mixcr assemble"
# mixcr assemble
mixcr assemble --write-alignments -f "$sample".IG.rescued2.vdjca "$sample".IG.clna

echo "mixcr contig assembly"
mixcr assembleContigs "$sample".IG.clna "$sample".IG.clns

echo "mixcr exportClones"
# mixcr exportClones
mixcr exportClones -count -fraction -vGene -dGene -jGene -vAlignment -dAlignment -jAlignment -aaFeature CDR3 -p fullImputed "$sample".IG.clns "$sample".IG.clones.tsv
# mixcr exportClones -f -o -t "$sample".IG.clns "$sample".IG.clones.tsv

sleep 5
echo "mixcr done"

echo "parse output started"

# parse output
cat "$sample".IG.clones.tsv | grep -v 'TRA' | grep -v 'TRB' | grep -v 'TRG' | grep -v 'TRD' > tmp.txt

# get all outputs
head -1 tmp.txt | sed -e "s/^/sample	/" > "$all_output"
tail -n +2 tmp.txt | sed -e "s/^/$id	/" >> "$all_output"

# get header
head -1 "$all_output" > "$sample".IG.clones.tsv.head1

# get top igh line
grep -m 1 "IGH" "$all_output" > "$sample".IG.clones.tsv.IGH

# get top igk line 
grep -m 1 "IGK" "$all_output" > "$sample".IG.clones.tsv.IGK

# get top igl line
grep -m 1 "IGL" "$all_output" > "$sample".IG.clones.tsv.IGL

# cat all files together
cat "$sample".IG.clones.tsv.head1 "$sample".IG.clones.tsv.IGH "$sample".IG.clones.tsv.IGK "$sample".IG.clones.tsv.IGL > "$top_output"

sleep 5
echo "parse output done"
