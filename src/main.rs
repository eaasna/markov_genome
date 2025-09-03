use clap::Parser;

mod args;
use crate::args::Args;

mod markov;
use crate::markov::run_markov_simulation;

fn main() {
    let mut args = Args::parse();

    if args.lens.len() == 0 {
        args.lens.push(1000);
    }
    
    run_markov_simulation(&args);
}

