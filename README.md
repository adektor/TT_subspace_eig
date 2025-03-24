# Inexact Subspace Projection Methods for Tensor Train Eigenvector Computation

This repository implements inexact subspace projection methods for computing approximate eigenvectors in the Tensor Train (TT) format.

## Requirements

This repository requires the [TT-Toolbox](https://github.com/oseledets/TT-Toolbox). Please ensure that it is installed and added to your MATLAB path before running any scripts.

## Usage

To compute approximate eigenvectors of a Heisenberg Hamiltonian using the Lanczos method and Chebyshev-filtered subspace iteration, run:

```matlab
compare_subspace
