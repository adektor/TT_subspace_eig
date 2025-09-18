# Inexact Subspace Projection Methods for Tensor Train Eigenvector Computation

Inexact Lanczos and subspace iteration methods for computing approximate eigenvectors in the Tensor Train (TT) format. 
https://arxiv.org/abs/2502.19578

## Requirements

This repository requires the [TT-Toolbox](https://github.com/oseledets/TT-Toolbox). Please ensure that it is installed and added to your MATLAB path before running any scripts. Run 

```matlab
setup
```
to add the necessary paths to your environment.

## Usage

To compare TT subspace iteration and TT Lanczos (reproduces results in Figure 6.1), run:

```matlab
subspace_krylov
```

To compare TT subspace iteration and DMRG (reproduces results in Figure 6.2), run:

```matlab
subspace_dmrg
```
