# data

## 5PSeq

Contains an input folder with the fasta files, gff files and meta data files required to run the paired_reads_pipeline.nf on the 5PSeq data set. Also, contains an placeholder output file to hold the ouputs of the pipeline.

## QuantSeq

Contains an input folder with the fasta files, gff files and meta data files required to run the paired_reads_pipeline.nf on the QuantSeq data set. Also, contains an placeholder output file to hold the ouputs of the pipeline.
This is for strains transformed with constructs:

- POT1-ccdB (empty vector control)
- pRPS3-mCherry-tRPS3_WT
- pRPS3-mCherry-tRPS3_mod0
- pRPS3-mCherry-tRPS3_modB
- pRPS3-mCherry-tRPS3_modE
- pTSA1-mCherry-tTSA1_WT
- pTSA1-mCherry-tTSA1_mod0
- pTSA1-mCherry-tTSA1_modB
- pTSA1-mCherry-tTSA1_modE
- pRPS3-mCherry-tRPS3_modA
- pTSA1-mCherry-tTSA1_modA
- pPGK1-mCherry-tRPS3_WT
- pPGK1-mCherry-tRPS3_mod0
- pPGK1-mCherry-tRPS3_modA
- pPGK1-mCherry-tRPS3_modB
- pPGK1-mCherry-tRPS3_modE
- pRPS3-mCherry-tRPS3_mod0


## QuantSeq_batch2

The same as above, but on the QuantSeq batch 2 dataset that included only constructs: 

- pRPS3-mCherry-tTOS6
- pRPS3-mCherry-tHSP26
- pRPS13-mCherry-tHSP26
- pRPS13-mCherry-tHSP26
- pHSP26-mCherry-tHSP26
- pRPS3-mCherry-tHSP26
- pRPS3-mCherry-tPMA1-long

## shared_data

Contains the fasta files and gff files for the raw plasmids and genome usead to create the reference files for each sample in both the 5PSeq and QuantSeq input folders. (Required to run `code/pre_processing/generate_genome_gff_and_fasta.Rmd``)
