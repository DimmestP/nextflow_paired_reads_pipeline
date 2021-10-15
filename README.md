# Template Nextflow Analysis Pipeline for the Analysis of Paired End Datasets
This GitHub repo holds the [NextFlow](https://www.nextflow.io/) workflow for reproducible analysis of the RNA-seq datasets created in the [Chimera project](https://github.com/DimmestP/chimera_project_manuscript).

We conducted two RNA-Seq assays to investigate alternative polyA site usage in the constructs created for the project:  

[QuantSeq](https://www.nature.com/articles/nmeth.f.376) is an RNA-Seq assay that quantifies the 3'UTR isoforms of mRNA transcripts. 

[5PSeq](https://www.nature.com/articles/nprot.2016.026) is an RNA-Seq assay that captures 5'phosphorylated mRNA transcript, enabling us to see the PolyA site usage of mRNA's flagged for decay. 

Each analysis has its own GitHub repo holding the specific data required for that assay but both a clones of this template repo to enable synchronisation of changes to core code.

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
