rm -r ${0}/construct_polyA_bam

mkdir ${0}/construct_polyA_bam

for f in ${0}/sorted_bam/*/*.bam;
do
	echo $f;
	sample_name=$(expr match $f '.*\([ABE][0-9]\+\)_')
        bedtools bamtobed -i $f |\
	bedtools intersect -a stdin -b\
	../data/${1}/input/genome_annotations/${sample_name}_sample_longest_full_ORF_with_constructs.gff | \
        grep -w - | grep  'mCherry\|YNL178W\|YML028W\|YCR012W' > ${0}/construct_polyA_bam/${sample_name}/${sample_name}_polyA_reads.bed
done;
