# code

### paired_reads_pipeline.nf

Nextflow script to filter, align and count RNA-seq reads. It is a template used in other analyses and requires several file name/location variables to be filled to run. 

```
# how to run code (for QuantSeq)
nextflow run paired_reads_pipeline.nf

# how to run code (for QuantSeq_heat_shock)
nextflow run paired_reads_pipeline.nf --combine_reads_over_lanes false --input_fq_dir '../data/QuantSeq_heat_shock/input/experiment_fastq/' --output_dir '../data/QuantSeq_heat_shock/output/' --sample_name_regex 'WD_[A-Z]\d+' --experiment_name 'QuantSeq_heat_shock'

# how to run code (for 5PSeq)
nextflow run paired_reads_pipeline.nf --experiment_name '5PSeq' --combine_reads_over_lanes false\
--input_fq_dir '/homes/wallacelab/datastore/wallace_rna/bigdata/fastq/5PSeq_data/'\
--output_dir '/homes/wallacelab/datastore/wallace_rna/data/2021/10-Oct/Sam/5PSeq_pipeline_output/'\
--fastq_file_regex '*dt_*.gz' --sample_name_regex '\d+[a-z]+' --read_1_forward true --remove_UMI true

```

### pre-processing/

Contains bash and R scripts used to format the input fasta, gff files, and index fasta files.

### post-processing/

Contains a bash script to extract a subset of reads mapped to genes of interest for further downstream analysis.
