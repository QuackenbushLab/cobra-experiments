# COBRA: Higher-order correction of persistent batch effects in correlation networks

Given a gene (co-)expression matrix and a design matrix, COBRA returns a decomposition of the covariance as a linear combination of components, one for each variable in the design matrix. 

COBRA can be applied for batch correction, covariate-specific co-expression analysis, and to study the effect of various parameters on the observed patterns of co-expression. 


<img align="center" width="100%" src="cobra.png">

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
- Marieke Kuijjer for the assistance in the ENCODE data pipeline. 
- Anahita Babaie, Rebekka Burkholz, Chen Chen, Derrick DeConti, Dawn DeMeo, Tara Eicher, Viola Fanfani, Jonas Fischer, Intekhab Hossain, Camila Lopes-Ramos, Panagiotis Mandros, John Quackenbush, Enakshi Saha, Katherine Shutta, and Ava Wilson for thoughtful critiques and discussions.

## Citation
If you find COBRA useful, star this repository!
