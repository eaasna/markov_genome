#!/bin/bash

set -ex

cd ../data

for k in 7 9 11 13 15; do
	cd k$k

	find . -name "markov*.fasta.shared.dump" | sort | xargs wc -l | awk '{print $1 "\t" $2}' | head -n -1  > tmp.counts

	find . -name "markov*.fasta" | sort | xargs head -n 1 | grep ">0_" | sed 's/>0_e//g' > tmp.errors
paste tmp.errors tmp.counts > markov.shared.counts 


	find . -name "mason*.fasta.shared.dump" | sort | xargs wc -l | awk '{print $1 "\t" $2}' | head -n -1  > tmp.counts

	find . -name "mason*.fasta" | sort | xargs head -n 1 | grep ">1_" | sed 's/>1_e//g' > tmp.errors
	paste tmp.errors tmp.counts > mason.shared.counts 
	
	cd ../
done
