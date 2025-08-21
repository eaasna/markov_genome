# Markov genome
Model a sequence database as a Markov process of order n and simulate random sequences with a similar k-mer count distribution. Markov genome ignores case and can model arbitrary 8-bit alphabets.  

## Download and Installation

<summary>Prerequisites </summary>

* Rust >= 1.74
* git

Refer to the [Rust Setup Tutorial](https://www.rust-lang.org/tools/install) for more in depth information.

<summary>Download current master branch</summary>

```bash
git clone git@github.com:eaasna/markov_genome.git
```
<summary>Building</summary>

```bash
cd markov_genome
cargo build --release
```
The `markov_genome` binary can be found in `target/release/`.

You may want to add the executable to your PATH:
```
export PATH=$(pwd)/target/release/markov_genome:$PATH
markov_genome --help
```

## Sequence simulation

Learn the transition probabilities of the test reference database and use a Markov process of order 3 to simulate two chromosomes of length 100bp:   
```
markov_genome --input test/ref.fasta --output sim.fasta --order 3 --lens 100 --lens 100
```

For a detailed list of options, see the help page:
```
markov_genome --help
```
