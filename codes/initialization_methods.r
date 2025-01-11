init_w_USVT = function(A, M, k = 3, tau = 1.01, upper = 0.5, lower = 5e-4){
	# -----------------------------------------------------------------------------
    # Description: Calculate initial estimation of Theta and P with hyperbolic model. 
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # A         - Adjacency matrix
    # M         - Matrix of masked entries
    # k         - Dimension of latent space
    # tau       - Thresholding parameter
    # upper     - upper bound for thresholding parameter
    # lower     - lower bound for thresholding parameter
    # -----------------------------------------------------------------------------
    # Output: 
    # Theta_0   - Initial estimation of Theta matrix
    # P_0       - Initial estimation of P matrix
    # -----------------------------------------------------------------------------

	svd_res = svd(A * M)
	n = dim(A)[1]
	if (is.null(k)){
		k = sum(svd_res$d > tau * sqrt(n * mean(A)))
	}
	# estimate P_0
	P_0 = svd_res$u[, 1 : k] %*% diag(svd_res$d[1 : k]) %*% t(svd_res$v[, 1 : k]) / mean(M)
	P_0[P_0 <= lower] = lower
	P_0[P_0 > upper] = upper
	P_0 = (P_0 + t(P_0)) / 2
	diag(P_0) = rep(0, n)

	# derive Theta_0 = - D_0
	Theta_0 = qlogis(P_0 / 2)
	diag(Theta_0) = rep(0, n)
	return(list(Theta_0 = Theta_0, P_0 = P_0))
}

###################################################################

init_KZ = function(Theta_0, K = 1, d = 3, Rmax = 4){
	# -----------------------------------------------------------------------------
    # Description: Calculate initial estimation of Z with hyperbolic model. 
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # Theta_0   - Estimation of Theta matrix
    # K         - Curvature parameter
    # d         - Dimension of latent space
    # Rmax      - The maximum distance of projected point to disk center
    # -----------------------------------------------------------------------------
    # Output: 
    # Z_0       - Initial estimation of latent positions
    # -----------------------------------------------------------------------------

	# G_{ij} = Z_i L Z_j^T
	G_0 = - cosh(- sqrt(K) * Theta_0)
	# bound entries in G_0 such that points are in a disk of given radius 
	radius = 3
	lower_bound = - cosh(2 * radius * sqrt(K))
	G_0[G_0 <= lower_bound] = lower_bound

	# decompose G to get Z following paper: hyperbolic distance matrix
	eigen_G = eigen(G_0)
	Lambda = diag(sqrt(abs(c(eigen_G$values[1 : (d - 1)], eigen_G$values[length(eigen_G$values)]))))
	Z_0 = t(diag(c(rep(1, d - 1), -1)) %*% Lambda %*% t(eigen_G$vectors[, c(1 : (d - 1), length(eigen_G$values))]))
	Z_0 = project_E_to_H(Z_0, d = d, Rmax = cosh(sqrt(K) * Rmax))
	return(Z_0)
}

###################################################################

init_EU = function(Theta_0, k = 2){
	# -----------------------------------------------------------------------------
    # Description: Calculate initial estimation of Z with Euclidean model. 
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # Theta_0   - Estimation of Theta matrix
    # k         - Dimension of latent space
    # -----------------------------------------------------------------------------
    # Output: 
    # Z_0       - Initial estimation of latent positions
    # -----------------------------------------------------------------------------
	
	n = dim(Theta_0)[1]
	D2_0 = Theta_0 * Theta_0

	# bound entries in D2_0 such that points are in a disk of given radius
	radius = 3
	D2_0[D2_0 >= radius ^ 2] = radius^2
	one_diag = rep(1, n) %*% t(D2_0[, 1])
	G_0 = - 0.5 * (D2_0 - one_diag - t(one_diag))
	eigen_G_0 = eigen(G_0)
	L_0 = diag(c(sort(sqrt(abs(eigen_G_0$values)), decreasing = TRUE)[1 : k]))
	Z_0 = eigen_G_0$vectors[, 1: k] %*% L_0
	return(Z_0)
}



