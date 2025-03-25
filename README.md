# Inexact Subspace Projection Methods for Tensor Train Eigenvector Computation

Inexact Lanczos and subspace iteration methods for computing approximate eigenvectors in the Tensor Train (TT) format. 
https://arxiv.org/abs/2502.19578

## Requirements

This repository requires the [TT-Toolbox](https://github.com/oseledets/TT-Toolbox). Please ensure that it is installed and added to your MATLAB path before running any scripts.

## Usage

To compute approximate eigenvectors of a Heisenberg Hamiltonian using the Lanczos method and Chebyshev-filtered subspace iteration and compare convergence, run:

```matlab
compare_subspace
