#!/bin/bash

set -ex 

markov_genome="../target/release/markov_genome"
k=9
kmc_dir="/home/evelin/bin"

run () {
	for k in 7 9 11 13 15; do
		mkdir -p ../data/k$k
		$kmc_dir/kmc -fm -k$k -ci1 -cs10000 ../data/$1 ../data/k$k/$1.res ../data
		for e in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8; do
			for seed in 11 22 33 44; do
				sim="../data/k$k/$2_e${e}_s${seed}.fasta"
				$markov_genome mutate --input ../data/$1 -e $e --seed $seed --output $sim
				$kmc_dir/kmc -fm -k$k -ci1 -cs100000 $sim $sim.res ../data/k$k
	        		$kmc_dir/kmc_tools simple ../data/k$k/$ref.res -ci1 -cx100000 $sim.res -ci1 -cx100000 intersect $sim.shared.res
				$kmc_dir/kmc_dump $sim.shared.res $sim.shared.dump
				rm $sim.res.kmc_* $sim.shared.res.kmc_* 
			done
		done
	done
	}

ref="markov_o3_fly.fasta"
prefix="markov"
run $ref $prefix

ref="mason.fasta"
prefix="mason"
run $ref $prefix

ref="query_concat.fasta"
prefix="fly"
run $ref $prefix
