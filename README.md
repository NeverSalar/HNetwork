HNetwork
================

This GitHub repository implements the estimation and inference methods in the paper "Hyperbolic Network Latent Space Model
with Learnable Curvature". 

## Getting Started

### Prerequisites

The implementation of the proposed method is built on top of the following packages in R (version 4.3.2).

* Rcpp 1.0.12
* RcppArmadillo 0.12.8.4.0
* RcppDist 0.1.1
* uniformly 0.5.0

## Usage
* HNetwork_1.0.tar.gz: The R package that implements the manifold gradient descent optimization algorithm, need to be installed before running experiments.
  
* simulation_script.r: Include pipeline codes for generating data, model estimation, and model evaluation.
* generate_data.r: Include functions for generating data from hyperbolic models and Euclidean models.
* initialization_methods.r: Include functions for calculating initial estimations of latent space curvature and positions.
* hyperbolic_toolbox.r: Include helper functions for visualization, geometric transformation, and sampleing in hyperbolic space.
