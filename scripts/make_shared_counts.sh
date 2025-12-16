#!/bin/bash

set -ex

data_dir="/buffer/ag_abi/evelina/simulated"
cd $data_dir


run () {
	for k in 5 7 9 11;do
		cd k$k

		name="$1"
		chr="$2"

		wc -l "$name.fa.dump" | awk -v n="$name" '{print $1 "\t./"n"_e0_s0.fasta.orig.dump" }' | head -n 1 > tmp.counts
		find . -name "$name*.fasta.shared.dump" | sort | xargs wc -l | awk '{print $1 "\t" $2}' | head -n -1  >> tmp.counts
		
		echo "0" > tmp.errors
		find . -name "$name*.fasta" | sort | xargs head -n 1 | grep '^>' | sed "s/>${chr}_e//g" >> tmp.errors
		paste tmp.errors tmp.counts > $name.shared.counts 
		
		cd ../
	done
}

run "markov" "0"
run "mason" "1"

run "mouse_chr4" "chr4"	
run "human_chr4" "chr4"

run "mouse_chr14" "chr14"	
run "human_chr14" "chr14"

if [ 0 -eq 1 ]; then

	for k in 5 7 9 11; do
	cd k$k
	find . -name "markov*.fasta.shared.dump" | sort | xargs wc -l | awk '{print $1 "\t" $2}' | head -n -1  > tmp.counts
	find . -name "markov*.fasta" | sort | xargs head -n 1 | grep ">0_" | sed 's/>0_e//g' > tmp.errors
paste tmp.errors tmp.counts > markov.shared.counts 


	find . -name "mason*.fasta.shared.dump" | sort | xargs wc -l | awk '{print $1 "\t" $2}' | head -n -1  > tmp.counts
	find . -name "mason*.fasta" | sort | xargs head -n 1 | grep ">1_" | sed 's/>1_e//g' > tmp.errors
	paste tmp.errors tmp.counts > mason.shared.counts 
	
	find . -name "mouse_chr4*.fasta.shared.dump" | sort | xargs wc -l | awk '{print $1 "\t" $2}' | head -n -1  > tmp.counts
	find . -name "mouse_chr4*.fasta" | sort | xargs head -n 1 | grep ">chr4_partial" | sed 's}>chr4_partial_e}}g' > tmp.errors
	paste tmp.errors tmp.counts > mouse_chr4.shared.counts 
	
	find . -name "human_chr4*.fasta.shared.dump" | sort | xargs wc -l | awk '{print $1 "\t" $2}' | head -n -1  > tmp.counts
	find . -name "human_chr4*.fasta" | sort | xargs head -n 1 | grep ">chr4_partial" | sed 's}>chr4_partial_e}}g' > tmp.errors	
	paste tmp.errors tmp.counts > human_chr4.shared.counts
	

	find . -name "mouse_chr14*.fasta.shared.dump" | sort | xargs wc -l | awk '{print $1 "\t" $2}' | head -n -1  > tmp.counts
	find . -name "mouse_chr14*.fasta" | sort | xargs head -n 1 | grep ">chr14_partial" | sed 's}>chr14_partial_e}}g' > tmp.errors
	paste tmp.errors tmp.counts > mouse_chr14.shared.counts 
	
	find . -name "human_chr14*.fasta.shared.dump" | sort | xargs wc -l | awk '{print $1 "\t" $2}' | head -n -1  > tmp.counts
	find . -name "human_chr14*.fasta" | sort | xargs head -n 1 | grep ">chr14_partial" | sed 's}>chr14_partial_e}}g' > tmp.errors	
	paste tmp.errors tmp.counts > human_chr14.shared.counts

	#find . -name "fly*.fasta.shared.dump" | sort | xargs wc -l | awk '{print $1 "\t" $2}' | head -n -1  > tmp.counts

	#find . -name "fly*.fasta" | sort | xargs head -n 1 | grep ">NC_004354.4" | sed 's}>NC_004354.4_e}}g' > tmp.errors
	#paste tmp.errors tmp.counts > fly.shared.counts 
		
	cd ../
done
fi

