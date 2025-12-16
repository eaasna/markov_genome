//! # Markov Genome
//!
//! `markov_genome` is a collection of utilities for applied Markov modelling in the context of 
//! biological sequence simulation. 

/// # Examples
/// markov_genome --help
/// markov_genome simulate --input test/ref.fasta --output sim.fasta --lens 100 --lens 200
/// markov_genome mutate --input sim.fasta --output mut.fasta --error 0.2 

use clap::Parser;

mod args;
use crate::args::{Cli, Commands};

mod markov;
use crate::markov::run_markov_simulation;

mod mutation;
use crate::mutation::run_mutation;

mod io;

fn main() {
    let args = Cli::parse();
 
    match args.command {
        Commands::Simulate (mut sim_args) => {
            println!("Simulating");
            if sim_args.lens.len() == 0 {
                sim_args.lens.push(1000);
            }
            run_markov_simulation(&sim_args);
        }
        Commands::Mutate (mut_args) => {
            println!("Mutating");
            run_mutation(&mut_args);
        }
    }
}
