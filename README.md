# lkv -- Linearized Kashiwara-Vergne
Code for the numerical computation of the dimension of the linearized Kashiwara-Vergne Lie algebra in various depths, see [[arXiv:2508.0808]](https://arxiv.org/abs/2508.08081).


# Prerequisites

You need a running rust installation, and Wolfram Mathematica.
Then just clone or download the repository.


# Usage
The computation of the dimension consists of two parts as explained in the paper: The computation of a lower bound and the computation of an upper bound.

## Upper bound 

For the upper bound only the rust code is required. Example usage from the command line, assuming that you are in the repository's base directory:

```sh
cargo run --release -- 12 20 -t 5 -d 1 -D 10
```

This computes (upper bounds on the) the dimensions of the linearized kv algebra in weights ("length") 12, 13, 14, 15, and depth from 1 to 10, using up to 5 threads.
The last three arguments are optional.

## Lower bound

The computation of the lower bound uses two steps. First, we need to compute formulas for the Brown polynomials. This is done using the mathematica notebook "brown_poly.nb". Execute the code in there. This should create files in the form "gXX_Y.json" that contain the Brown generators.
If the files are present, then run the rust code, for example:

```sh
cargo run --release -- 12 20 --lower -t 5 -d 1 -D 10
```

This computes lower bounds for the dimension of the extended (!) linearized KV algebra. 
Mind that only non-zero entries are reported in the summary displayed by the program.
Also mind that in order to compute the lower bound up to weight W the mathematica code has to be run up to that weight beforehand.

## On the ranges in the paper

Generally, to obtain the full range of entries listed in the paper significant computational resources are necessary.
Also, depending on the weight and depth range, we used different integer datatypes to save memory.
Particularly, for low weights but high depths the type CycInt in cyc_and_lie.rs should be changed to smaller integers to conserve memory.
For example, for weights up to 31 it is sufficient to use u32 etc.

