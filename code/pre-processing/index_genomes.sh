#!/bin/bash

for file in ../../data/${1}/input/genome_fasta/*.fa
do
	hisat2-build $file "../../data/${1}/input/indexed_genome/$(basename -s .fa $file)" 
done
