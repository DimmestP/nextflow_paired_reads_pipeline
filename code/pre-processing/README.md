# preprocessing

### generate_genome_gff_and_fasta.Rmd
R markdown file combining reference genome sequence and annotation with construct specific plasmid sequences and annotations. 

### index_genomes.sh
Bash script for indexing all fasta files in data/input/genome_fasta/pipeline_fasta. Requires HiSat2 installed to run

```
# how to run script
# enable script to be executed
chmod +x index_genomes.sh

# run script
./index_genomes.sh
```
