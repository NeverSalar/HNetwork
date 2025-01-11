library(Rcpp)
library(RcppArmadillo)
library(RcppDist)
library(uniformly)
library(HNetwork)

source("codes/hyperbolic_toolbox.r")
source("codes/generate_data.r")
source("codes/simulation_complex.r")
source("codes/initialization_methods.r")

sim = function(n = 500, k = 3, K = 1, R = pi, init_K, eta_K = 5 / n / n, eta_Z = 5 / n, eps = 1e-6, eps_K = 1e-8, max_iter = 1000, seed = NULL){
	# -----------------------------------------------------------------------------
    # Description: Example of a simulation pipeline, from data generation to estimation to evaluation.
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # n         - Number of nodes
    # k         - Dimension of latent space
    # K         - Curvature parameter
    # R         - Radius
    # init_K    - Initial value of curvature parameter
    # eta_K     - Stepsize for updating K
    # eta_Z     - Stepsize for updating Z
    # eps       - Threshold parameter for loss
    # eps_K     - Threshold parameter for K
    # max_iter  - Maximum number of iteration 
    # seed      - Random seed for generating data

    # -----------------------------------------------------------------------------
    # Output: 
    # esti_H2, esti_H3, esti_E2  - Estimations, including estimation of latent positions and curvature
    # metrics_H2, metrics_H3, metrics_E2  - Estimation errors for K, Z, P, Theta
    # -----------------------------------------------------------------------------
	
    # entries observable by algorithm, can be masked for link prediction experiments
    M = matrix(1, n, n)

    dat = gen_data(n = n, K = K, R = R, seed = seed)
    # initialize with higher dimension Euclidean model
    init_TP = init_w_USVT(dat$A, k, 10)
    init_Z = init_EU(init_TP$Theta_0, 20)
    esti_EU20 = GD_cpp_EU(dat$A, M, 20, dat$P, dat$Theta, init_Z, eta_Z, eps, max_iter)

    # H2
    init_Z = init_KZ(esti_EU20$Theta, K = init_K, d = 3, 10)
    esti_H2 = GD_cpp(dat$A, M, 3, dat$Z[, c(1,2,3)], dat$P, dat$Theta, dat$K, init_K, init_Z, eta_K, eta_Z, eps, max_iter)
    metrics_H2 = c(esti_H2$ERR_P[esti_H2$iter], esti_H2$ERR_Theta[esti_H2$iter], esti_H2$ERR_MZ[esti_H2$iter], esti_H2$ERR_K[esti_H2$iter])

    # H3
    init_Z = init_KZ(esti_EU20$Theta, K = init_K, d = 4, 10)
    esti_H3 = GD_cpp(dat$A, M, 4, dat$Z, dat$P, dat$Theta, dat$K, init_K, init_Z, eta_K, eta_Z, eps, max_iter)
    metrics_H3 = c(esti_H3$ERR_P[esti_H3$iter], esti_H3$ERR_Theta[esti_H3$iter], esti_H3$ERR_MZ[esti_H3$iter], esti_H3$ERR_K[esti_H3$iter])

    # E2
    init_Z = init_EU(init_TP$Theta_0, 2)
    esti_E2 = GD_cpp_EU(dat$A, M, 2, dat$P, dat$Theta, init_Z, eta_Z, eps, max_iter)
    metrics_E2 = c(esti_E2$ERR_P[esti_E2$iter], esti_E2$ERR_Theta[esti_E2$iter], esti_E2$ERR_MZ[esti_E2$iter], esti_E2$ERR_K[esti_E2$iter])

    res = list(esti_H2 = esti_H2, esti_H3 = esti_H3, esti_E2 = esti_E2, metrics_H2 = metrics_H2, metrics_H3 = metrics_H3, metrics_E2 = metrics_E2)
    return(res)
    }






