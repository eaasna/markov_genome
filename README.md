# Markov genome
Model a sequence database as a Markov process of order n and simulate random sequences with a similar k-mer count distribution.

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
