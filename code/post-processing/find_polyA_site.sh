rm -r ${1}/construct_polyA_bam

mkdir ${1}/construct_polyA_bam

# create regex for genes of interest
regex='mCherry\|YNL178W\|YML028W\|YCR012W'

# create lamda function to filter gff files for only entries relating to genes of interest
gff_genes_of_interest="grep $regex ../../data/${2}/input/genome_annotations/${sample_name}_sample_longest_full_ORF_with_constructs.gff"

for f in ${1}/sorted_bam/*/*.bam;
do
	echo $f;
	sample_name=$(expr match $f '.*\([ABE][0-9]\+\)_')
        mkdir ${1}/construct_polyA_bam/${sample_name}/
        bedtools bamtobed -i $f |\
	bedtools intersect -u -a stdin -b $gff_genes_of_interest | \
        grep -w - > ${1}/construct_polyA_bam/${sample_name}/${sample_name}_polyA_reads.bed
done;
