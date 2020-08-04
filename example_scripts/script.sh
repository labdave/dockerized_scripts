# this is a simple script that counts the number of lines in a bam header

bam=""
while getopts ":b::o::" opt; do
  case $opt in
    b) bam="$OPTARG"
    ;;
    o) output="$OPTARG"
	;;
    \?) echo "Invalid option -$OPTARG" >&2
    ;;
  esac
done

if [ "$bam" != "" ] && [ "$output" != "" ]
then
	samtools view -H $bam | wc -l > $output
else
	echo "Hello from inside script.sh"
fi