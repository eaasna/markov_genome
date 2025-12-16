#!/bin/bash

set -ex 

markov_genome="../target/release/markov_genome"
kmc_dir="/home/evelin/bin"
kmc_dir="/group/ag_abi/evelina/decode/bin"
data_dir="/buffer/ag_abi/evelina/simulated"
run () {
	#for k in 7 9 11 13 15 17; do
	for k in 5 7 9 11; do
		mkdir -p $data_dir/k$k
		$kmc_dir/kmc -fm -k$k -ci1 -cs10000 $data_dir/$1 $data_dir/k$k/$1.res $data_dir
		$kmc_dir/kmc_dump $data_dir/k$k/$1.res $data_dir/k$k/$1.dump

		for e in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8; do
			for seed in 11 22 33 44; do
				sim="$data_dir/k$k/$2_e${e}_s${seed}.fasta"
				$markov_genome mutate --input $data_dir/$1 -e $e --seed $seed --output $sim
				
				$kmc_dir/kmc -fm -k$k -ci1 -cs100000 $sim $sim.res $data_dir/k$k
				$kmc_dir/kmc_tools simple $data_dir/k$k/$ref.res -ci1 -cx100000 $sim.res -ci1 -cx100000 intersect $sim.shared.res
				$kmc_dir/kmc_dump $sim.shared.res $sim.shared.dump
				rm $sim.res.kmc_* $sim.shared.res.kmc_* 
			done
		done
	done
	}

ref="markov.fa"
prefix="markov"
run $ref $prefix

ref="mason.fa"
prefix="mason"
run $ref $prefix

ref="mouse_chr4.fa"
prefix="mouse_chr4"
run $ref $prefix

ref="human_chr4.fa"
prefix="human_chr4"
run $ref $prefix

ref="mouse_chr14.fa"
prefix="mouse_chr14"
run $ref $prefix

ref="human_chr14.fa"
prefix="human_chr14"
run $ref $prefix

############################################################################################

#ref="markov_order3_fly.fa"
#prefix="markov"
#run $ref $prefix

#ref="mason_fly.fa"
#prefix="mason"
#run $ref $prefix

#ref="fly.fa"
#prefix="fly"
#run $ref $prefix

#ref="mouse_140Mb.fa"
#prefix="mouse_chr4"
#run $ref $prefix

#ref="human_140Mb.fa"
#prefix="human_chr7"
#run $ref $prefix
