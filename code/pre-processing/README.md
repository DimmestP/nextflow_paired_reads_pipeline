# pre-processing

### generate_genome_gff_and_fasta.Rmd
R markdown file combining reference genome sequence and annotation with construct specific plasmid sequences and annotations. Creates QuantSeq references by default, change experiment_name to 5PSeq in first code block if needed.

```
# how to run script
# open R in console
R

# run Rmd with knitr
knitr::knit("generate_genome_gff_and_fasta.Rmd")
``` 

### index_genomes.sh
Bash script for indexing all fasta files in data/input/genome_fasta/pipeline_fasta. Requires HiSat2 installed to run

```
# how to run script
# enable script to be executed
chmod +x index_genomes.sh

# run script passing whether the experiment is QuantSeq or 5PSeq 
# so that it checks correct file of .fa files 
./index_genomes.sh QuantSeq
```
