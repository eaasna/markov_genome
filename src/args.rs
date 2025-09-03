use clap::Parser;

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
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

