# Markov genome [![install with bioconda][1]][2]
Model a sequence database as a Markov process of order n and simulate random sequences with a similar k-mer count distribution. Markov genome ignores case and can model arbitrary 8-bit alphabets.  

[1]: https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat "Install with bioconda"
[2]: #install-with-bioconda-linux
<!--
    This is the CI badge image:
        `https://img.shields.io/github/workflow/status/` - we do not use GitHub's badges as they are not customisable.
        `/seqan/dream-stellar/` - owner/repository
        `CI%20on%20Linux` - name of the workflow as encoded URL (e.g., whitespace = %20)
        `main` - branch to show
        `?style=flat&logo=github` - use a GitHub-style badge
        `&label=markov_genome%20CI` - text on the badge
        `"Open GitHub actions page"` - this text will be shown on hover
-->


## Download and Installation
### Install with [cargo](https://crates.io/crates/markov_genome) (Linux)

```bash
cargo install markov_genome
```

### Install with [bioconda](http://bioconda.github.io/recipes/markov_genome/README.html) (Linux)

```bash
conda install -c bioconda -c conda-forge markov_genome
```
### Build from source

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
