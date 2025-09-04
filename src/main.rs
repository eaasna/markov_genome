use clap::Parser;
use clap::{arg, command, Command, value_parser};
use std::path::PathBuf;

mod args;
use crate::args::{Cli, Commands};

mod markov;
use crate::markov::run_markov_simulation;

mod mutation;
use crate::mutation::run_mutation;

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
