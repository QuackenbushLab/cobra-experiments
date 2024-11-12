# COBRA: Higher-order correction of persistent batch effects in correlation networks
**Published as [original paper](https://academic.oup.com/bioinformatics/article/40/9/btae531/7748404) on Bioinformatics, Volume 40, Issue 9, September 2024**

## Abstract

Systems biology analyses often use correlations in gene expression profiles to infer co-expression networks that are then used as input for gene regulatory network inference or to identify functional modules of co-expressed or putatively co-regulated genes. While systematic biases, including batch effects, are known to induce spurious associations and confound differential gene expression analyses (DE), the impact of batch effects on gene co-expression has not been fully explored. Methods have been developed to adjust expression values, ensuring conditional independence of mean and variance from batch or other covariates for each gene, resulting in improved fidelity of DE analysis. However, such adjustments do not address the potential for spurious differential co-expression (DC) between groups. Consequently, uncorrected, artifactual DC can skew the correlation structure, leading to the identification of false, non-biological associations, even when the input data are corrected using standard batch correction.

In this work, we demonstrate the persistence of confounders in covariance after standard batch correction using synthetic and real-world gene expression data examples. We then introduce Co-expression Batch Reduction Adjustment (COBRA), a method for computing a batch-corrected gene co-expression matrix based on estimating a conditional covariance matrix. COBRA estimates a reduced set of parameters expressing the co-expression matrix as a function of the sample covariates, allowing control for continuous and categorical covariates. COBRA is computationally efficient, leveraging the inherently modular structure of genomic data to estimate accurate gene regulatory associations and facilitate functional analysis for high-dimensional genomic data.

## Overview

<img align="center" width="90%" src="cobra.png">

Given a gene (co-)expression matrix and a design matrix, COBRA returns a decomposition of the covariance as a linear combination of components, one for each variable in the design matrix. 

COBRA can be applied for batch correction, covariate-specific co-expression analysis, and to study the effect of various parameters on the observed patterns of co-expression. 

## Usage

In ```netbooks``` we provide a simple tutorial to show how to use COBRA in different applications. Check-out the [R Version](https://github.com/netZoo/netbooks/blob/main/netbooks/netZooR/COBRA.ipynb) and [Python Version](https://github.com/netZoo/netbooks/blob/main/netbooks/netZooPy/cobra.ipynb). For an interactive playground [click here](http://netbooks.networkmedicine.org/hub/login)!

COBRA is part of the [Network Zoo](https://netzoo.github.io/), and the source code is available both in and [netZooR](https://github.com/netZoo/netZooR) and [netZooPy](https://github.com/netZoo/netZooPy).

## Minimal examples

### Python

```python
from netZooPy.cobra import *
import numpy as np

# Generate toy data
G = 1000 # Genes
n = 100 # Samples (e.g. individuals)
q = 2 # Covariates in the design matrix

# expression of size (G, n); design matrix of size (G, q)
expression = np.random.random((G, n))
X = np.vstack(([1 for i in range(n)], np.random.rand(n))).T # The first column of X is an intercept

# Run COBRA
psi, Q, d, g = cobra(X, expression)
```

### R

```r
library(netZooR)

# Generate toy data
G <- 1000 # Genes
n <- 100 # Samples (e.g. individuals)
q <- 2 # Covariates in the design matrix

# expression of size (G, n); design matrix of size (G, q)
expression <- matrix(runif(G * n), nrow = G)
X <- cbind(rep(1, n), runif(n)) # The first column of X is an intercept

# Run COBRA
cobra_obj <- cobra(X, estimates)

# Extract matrix Q and covariate-specific eigenvalues psi
Q <- cobra_obj$Q
psi <- cobra_obj$psi
```

## Structure of the repo

- ```data``` contains the data used and generated in the experiments. Note that the version on GitHub does **not** contain all the data. Please download them from [Zenodo](https://zenodo.org/records/10154758). 
- ```experiments``` contains the source code to reproduce our experiments.
- ```figures``` are the output figures from the experiments.

## Appreciation
- Anahita Babaie for useful code optimizations.
- Marieke Kuijjer for the assistance in the ENCODE data pipeline.
- Anahita Babaie, Rebekka Burkholz, Chen Chen, Derrick DeConti, Dawn DeMeo, Tara Eicher, Viola Fanfani, Jonas Fischer, Intekhab Hossain, Camila Lopes-Ramos, Panagiotis Mandros, John Quackenbush, Enakshi Saha, Katherine Shutta, and Ava Wilson for thoughtful critiques and discussions.

## Citation
```bibtex
@article{micheletti2024higher,
  title={Higher-order correction of persistent batch effects in correlation networks},
  author={Micheletti, Soel and Schlauch, Daniel and Quackenbush, John and Ben Guebila, Marouen},
  journal={Bioinformatics},
  volume={40},
  number={9},
  pages={btae531},
  year={2024},
  publisher={Oxford University Press}
}
```

If you find COBRA useful, star this repository!
