use std::ffi::OsStr;
use std::ffi::OsString;
use std::path::PathBuf;

use clap::{Args, Parser, Subcommand, ValueEnum};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct SimulateArgs {
    #[arg(short, long)]
    pub input: String,

    #[arg(short, long, default_value_t = String::from("seq.fasta"))]
    pub output: String,
    
    #[arg(short, long)]
    pub lens: Vec<usize>,

    #[arg(long, default_value_t = 3)]
    pub order: usize,

    #[arg(long, default_value_t = 42)]
    pub seed: u64,

    #[arg(short, long, default_value_t = false)]
    pub verbose: bool,
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct MutateArgs {
    #[arg(short, long)]
    pub input: String,

    #[arg(short, long, default_value_t = String::from("seq.fasta"))]
    pub output: String,
    
    #[arg(long, default_value_t = 42)]
    pub seed: u64,

    #[arg(short, long, default_value_t = 0.1)]
    pub error: f64,

    #[arg(short, long, default_value_t = false)]
    pub verbose: bool,
}
#[derive(Debug, Parser)] // requires `derive` feature
#[command(name = "markov_genome")]
#[command(about = "Markov chain sequence simulation", long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Debug, Subcommand)]
pub enum Commands {
    #[command(arg_required_else_help = true)]
    Simulate (SimulateArgs),
    Mutate (MutateArgs),
}

