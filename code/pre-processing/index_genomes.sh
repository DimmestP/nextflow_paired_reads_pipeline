#!/bin/bash

for file in ../../data/QuantSeq/input/genome_fasta/*.fa
do
	hisat2-build $file "../../data/QuantSeq/input/indexed_genome/$(basename -s .fa $file)" 
done
