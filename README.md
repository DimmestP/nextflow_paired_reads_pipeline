# Template Nextflow Analysis Pipeline for the Analysis of Paired End Datasets
This GitHub repo holds the [NextFlow](https://www.nextflow.io/) workflow for reproducible analysis of the RNA-seq datasets created in the [Chimera project](https://github.com/DimmestP/chimera_project_manuscript).

We conducted two RNA-Seq assays to investigate alternative polyA site usage in the constructs created for the project:  

[QuantSeq](https://www.nature.com/articles/nmeth.f.376) is an RNA-Seq assay that quantifies the 3'UTR isoforms of mRNA transcripts. 

[5PSeq](https://www.nature.com/articles/nprot.2016.026) is an RNA-Seq assay that captures 5'phosphorylated mRNA transcript, enabling us to see the PolyA site usage of mRNA's flagged for decay. 

Each analysis can be conducted using the code in this repo if the appropriate user defined variables are passed to the nextflow pipeline (see code folder README).

# Installation
First, follow the installation steps on the [NextFlow](https://www.nextflow.io/).

To run this pipeline you will also need to have [conda](https://www.nextflow.io/docs/latest/index.html) installed as each process runs in a conda environment.

Then, clone this repo to your local computer/cluster. 

# Usage
### Main analysis
To run this pipeline move to the code folder and pass the paired_reads_pipeline.nf file to nextflow. By default, the pipeline will analyse the QuantSeq data held on the lab datastore, see the README file in the code folder to see this.

```
cd code
nextflow run paired_reads_pipeline.nf
```

### Pre-processing
The fasta files and gff files which which the pipeline will align and count reads to needs to be created before running the pipeline. The reference fasta files also need to be properly index before HISAT2 begin align. And Rmd file creates the required fasta/gff files and a bash script indexes the fasta files. To run these scripts for the quantseq pipeline see below (more details are availble in the code/pre-processsing folder README file for analysing 5PSeq).

```
# create fasta and gff files
# move to the correct folder
cd ./code/pre-processing

# start the R cli
R

# run the Rmd file
knitr::knit("generate_genome_gff_and_fasta.Rmd")

# index the ref genome fastas
# enable script to be executed
chmod +x index_genomes.sh

# run script passing whether the experiment is QuantSeq or 5PSeq 
# so that it checks correct file of .fa files 
./index_genomes.sh QuantSeq
```

# Post-processing
The output bam files are quite large for each sample so I created a script to extract only the reads mapped to genes of interest for each sample (for easier downstream analysis).
```
# create fasta and gff files
# move to the correct folder
cd ./code/post-processing

# how to run script
# enable script to be executed
chmod +x find_polyA_site.sh

# run script whilst passing the path to the folder with the
# find_polyA_site nextflow pipeline output and the 
# experiment name
./find_polyA_site.sh /homes/wallacelab/datastore/wallace_rna/data/2021/10-Oct/Sam/chimera_quantseq_pipeline_output QuantSeq
```

# Status
The pipeline has the following functionality (Any unticked funtionality is currently being developed)

- [x] combineLanesAcrossSamples
- [x] Fastqc
- [x] cutAdapters
- [x] alignHisat2
- [x] samViewSort
- [x] makeBedgraphs
- [x] renameBamSample
- [x] countAllmRNA
- [x] runMultiQC
