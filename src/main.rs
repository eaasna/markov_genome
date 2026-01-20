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

mod sequence_model;

mod simulation;
use crate::simulation::run_simulation;

mod mutation;
use crate::mutation::run_mutation;
use crate::sequence_model::{SequenceGrammarModel, MarkovModel};

mod io;

mod grammar;

fn main() {
    let args = Cli::parse();
 
    match args.command {
        Commands::Simulate (mut sim_args) => {
            println!("Simulating");
            if sim_args.lens.len() == 0 {
                sim_args.lens.push(1000);
            }
            if sim_args.annotation.is_empty() {
                let sequence_model = MarkovModel::new(&sim_args);
                run_simulation(&sequence_model, &sim_args);
            } else {
                let sequence_model = SequenceGrammarModel::new(&sim_args);
                run_simulation(&sequence_model, &sim_args);
            }
        }
        Commands::Mutate (mut_args) => {
            println!("Mutating");
            run_mutation(&mut_args);
        }
    }
}
