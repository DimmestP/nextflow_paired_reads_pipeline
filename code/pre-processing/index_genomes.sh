#!/bin/bash

for file in ../../data/${1}/input/genome_fasta/*.fa
do
	~/bin/hisat2-2.1.0/hisat2-build $file "../../data/${1}/input/indexed_genome/$(basename -s .fa $file)" 
done
