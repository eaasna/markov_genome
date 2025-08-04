use clap::Parser;
use bio::io::fasta::{Reader};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
pub struct Args {
    #[arg(short, long)]
    input: String,

    #[arg(short, long, default_value_t = String::from("seq.fasta"))]
    output: String,
    
    #[arg(short, long, default_value_t = 100)]
    seq_len: usize,
}

fn main() {
    let args = Args::parse();

    let reader = Reader::from_file(args.input);
    let records = reader.expect("fasta reader: got an io::Error or could not read_line()").records();

    println!("Input FASTA");
    let mut rec_in_vec : Vec<Vec<char>> = Vec::new(); 
    for result in records {
        let record = result.expect("Error during fasta record parsing");
        let mut c_vec : Vec<char> = Vec::new();
        println!("{}", record.id());
        for i_val in record.seq() {
            let c_val = match i_val {
                b'A' | b'a' => 'A',
                b'C' | b'c' => 'C',
                b'G' | b'g' => 'G',
                b'T' | b't' => 'T',
                _ => panic!("Unknown char\t{:b}", i_val),
            //TODO: read char from int code
            };
            print!("{}", c_val);
            c_vec.push(c_val);
        }
        println!();
        rec_in_vec.push(c_vec);
    }

    
    use markov::Chain;
    let mut chain = Chain::of_order(3);
    // train the chain on some vectors
    for i in 0..rec_in_vec.len() {
        chain.feed(rec_in_vec[i].clone());
    }

    println!("Output FASTA");
    // generate sequences
    let mut rec_out_vec : Vec<Vec<char>> = Vec::new();
    
    let gen = chain.iter_for(args.seq_len);
    let mut rec_out : Vec<char> = Vec::new();
    for c in gen {
        rec_out.push(c[0]);
        println!("{}", rec_out.len())
    }
    rec_out_vec.push(rec_out.clone());
    
    for i in 0..rec_out_vec.len() {
        for j in 0..rec_out_vec[i].len() {
            print!("{}", rec_out_vec[i][j]);
        }
        println!();
    }
}

