#!/bin/bash

for o in 3 5 7; do
	../target/debug/markov_genome simulate --input ref.fasta --verbose --kmer 20 --annotation genmap.txt --output sim_o$o.fasta --order $o --lens 100 --lens 100
done

