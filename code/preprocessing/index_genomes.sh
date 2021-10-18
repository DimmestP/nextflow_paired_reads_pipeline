#!/bin/bash

for file in ../../data/input/genome_fasta/pipeline_fasta/*.fa
do
	hisat2-build $file "../../data/input/indexed_genome/$(basename -s .fa $file)" 
done
