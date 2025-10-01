#!/bin/bash

set -ex

mason="/group/ag_abi/evelina/DREAM-stellar-benchmark/lib/raptor_data_simulation/build/src/mason2/src/mason2-build/bin/mason_genome"

order=3
out_dir="/buffer/ag_abi/evelina/simulated"
seed=13

name="markov"
../target/release/markov_genome simulate --seed $seed --order $order --input /srv/data/evelina/fly/dna4.fasta --lens 10000 --output $out_dir/${name}.fa

name="mason"
$mason -l 10000 -o $out_dir/$name.fa --seed $seed

name="mouse_chr4"
echo ">chr4_partial" > $out_dir/$name.fa
head -n 7708188 /srv/data/evelina/mouse/dna4.fasta | tail -n 125 >> $out_dir/$name.fa

name="mouse_chr14"
echo ">chr14_partial" > $out_dir/$name.fa
head -n 24097823 /srv/data/evelina/mouse/dna4.fasta | tail -n 125 >> $out_dir/$name.fa

name="human_chr4"
echo ">chr4_partial" > $out_dir/$name.fa
head -n 9618075 /srv/data/evelina/human/dna4.fasta | tail -n 125 >> $out_dir/$name.fa

name="human_chr14"
echo ">chr14_partial" > $out_dir/$name.fa
head -n 27992611 /srv/data/evelina/human/dna4.fasta | tail -n 125 >> $out_dir/$name.fa
