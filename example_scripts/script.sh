# this is a simple script that counts the number of lines in a bam header

bam=""
while getopts ":b::" opt; do
  case $opt in
    b) bam="$OPTARG"
    ;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

if [ "$bam" != "" ]
then
	samtools view -H $bam | wc -l
else
	echo "Hello from inside script.sh"
fi