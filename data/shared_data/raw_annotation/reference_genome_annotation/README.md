# reference_genome_annotation

### longest_full-ORF_ypd_plus_other_fixed_UTR_length_transcripts.gff
Genome wide annotation for the S288C reference genome with each gene annotated with a 'primary_transcript' consisting of an ORF, 3'UTR and 5'UTR (As stated by Pelechano et al 2013).

### URA3_deletion_position.bed
Bed file with one entry corresponding to the position of the URA3 deletion in the BY4741 genome. Required for masking using the bedtools tool called in the code/pre-processing/generate_genome_gff_and_fasta.Rmd script.
