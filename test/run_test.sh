#!/bin/bash

for o in 3 5 7; do
	k=20
	../target/debug/markov_genome simulate --input ref.fasta --verbose --kmer $k --annotation genmap.txt --output sim_o${o}_k${k}.fasta --order $o --lens 100 --lens 100 &> annotated.out 
	../target/debug/markov_genome simulate --input ref.fasta --verbose --output sim_o$o.fasta --order $o --lens 100 --lens 100 &> simple.out

	ann_sum=$(sha256sum annotated.out | awk '{print $1}')
	sim_sum=$(sha256sum simple.out | awk '{print $1}')
	if [ $ann_sum != $sim_sum ]; then
		echo "Not equal\t"
	fi
done

