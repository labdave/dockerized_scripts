sample=$1
input_file=$2
threads=$3
thresh=$4
bed=$5
good_output=$6
bad_output=$7

# keep the header and greater than 500 insert sizes
echo "extracting reads with insert sizes greater than $thresh"
samtools view -@ $threads -h $input_file | awk -v var="$thresh" 'substr($0,1,1)=="@" || ($9>var) || ($9<-var)' | samtools view -@ $threads -b > $sample.is_gt_$thresh.bam

# intersect the gt500 bam file with an exon bed only allowing reads that fully map to exons and create a bed output
echo "intersect bam with pad 10 exon file"
bedtools intersect -f 1 -wao -bed -abam $sample.is_gt_$thresh.bam -b $bed > $sample.bed

# remove reads that have 0 length intersections
echo "select reads that actually intersect"
awk '($17>0)' $sample.bed > $sample.nonzero.bed

# run custom python script that finds read pairs that have intersections in different exons of the same gene
echo "obtain actual list of bad reads"
python3 remove_long_inserts.py $sample.nonzero.bed $sample.bad_reads.txt

# remove these bad reads from the original bam file
echo "filter samreads based on bad reads list"
java -jar picard.jar FilterSamReads I=$input_file O=$good_output RLF=$sample.bad_reads.txt FILTER=excludeReadList

# keep these bad reads from the original bam file
echo "keep samreads based on bad reads list"
java -jar picard.jar FilterSamReads I=$input_file O=$bad_output RLF=$sample.bad_reads.txt FILTER=includeReadList
