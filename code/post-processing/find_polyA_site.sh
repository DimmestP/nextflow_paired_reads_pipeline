rm -r ${1}/construct_polyA_bam

mkdir ${1}/construct_polyA_bam

samtool_TMPFILE=$(mktemp) || exit 1

gff_TMPFILE=$(mktemp) || exit 1

# create regex for genes of interest
regex="mCherry\|YNL178W\|YML028W\|YCR012W"

for f in ${1}/sorted_bam/*/*.bam;
do
	echo $f;
	sample_name=$(expr match $f '.*/\([ABEdt0-9]\+\)_')

	# create temp gff file to filter gff files for only entries relating to genes of interest
	grep $regex ../../data/${2}/input/genome_annotations/${sample_name}_sample_longest_full_ORF_with_constructs.gff > $gff_TMPFILE

        mkdir ${1}/construct_polyA_bam/${sample_name}/
	
        samtools view -F 396  $f |\
	grep NH:i:1 \
	> $samtool_TMPFILE
	
	bedtools intersect -u -S -bed -abam samtool_TMPFILE -b $gff_TMPFILE \
        > ${1}/construct_polyA_bam/${sample_name}/${sample_name}_polyA_reads.bed
done;
