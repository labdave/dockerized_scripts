# this is a simple script that counts the number of lines in a bam header

bam=""
while getopts ":b::l::c::" opt; do
  case $opt in
    b) bam="$OPTARG"
    ;;
    l) lines="$OPTARG"
    ;;
    c) chars="$OPTARG"
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

if [ "$bam" != "" ] && [ "$lines" != "" ] && [ "$chars" != "" ]
then
	samtools view -H $bam | wc -l > $lines
  samtools view -H $bam | wc -m > $chars
else
	echo "Hello from inside script.sh"
fi