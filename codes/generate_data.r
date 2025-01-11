gen_Z_3d = function(n = 100, K = 1, R = pi / 2){
	# -----------------------------------------------------------------------------
    # Description: Uniform sample latent positions on 3-dim Poincare ball.
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # n         - Number of nodes
    # R         - Radius
    # K         - Curvature parameter
    # -----------------------------------------------------------------------------
    # Output: 
    # (Z1, Z2, Z3, Z4) - Coordinates of latent positions on hyperboloid model
    # -----------------------------------------------------------------------------

	# Step 1: Sample radial distance
	u <- runif(n, 0, 1)
	r <- acosh(1 + u * (cosh(sqrt(K) * R) - 1)) / sqrt(K)

	# Step 2: Sample angles uniformly
	theta <- runif(n, 0, 2 * pi)
	phi <- acos(runif(n, -1, 1))

	# Step 3: Convert to hyperboloid model coordinates
	Z4 <- cosh(r)
	Z3 <- sinh(r) * sin(phi) * cos(theta)
	Z2 <- sinh(r) * sin(phi) * sin(theta)
	Z1 <- sinh(r) * cos(phi)

	return(cbind(Z1 = Z1, Z2 = Z2, Z3 = Z3, Z4 = Z4))
}

###################################################################

gen_data_3d = function(n = 100, K = 1, R = pi / 2, Z = NULL, seed = NULL, seed_Z = NULL){
	# -----------------------------------------------------------------------------
    # Description: Generate data for simulation study with 3-dim hyperbolic model.
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # n         - Number of nodes
    # K         - Curvature parameter
    # R         - Radius
    # seed      - Random seed for generating edges
    # seed_Z    - Random seed for generating latent positions

    # -----------------------------------------------------------------------------
    # Output: 
    # Z         - Matrix of latent positions
    # P         - Matrix of linkage probability
    # D         - Matrix of pairwise distance
    # Theta     - Matrix of theta parameter 
    # A         - Adjacency matrix
    # n         - Number of nodes
    # K         - Curvature parameter
    # R         - Radius
    # -----------------------------------------------------------------------------

	if (is.null(Z)){
		set.seed(seed_Z)
		Z = gen_Z_3d(n = n, K = 1, R = R * sqrt(K))
	}
	DPT = gen_P(n = n, K = K, R = R, Z = Z, alpha = alpha)
	P = DPT$P
	D = DPT$D
	Theta = DPT$Theta
	set.seed(seed)
	A = gen_network(n = n, K = K, R = R, P = P, seed = seed)
	return(list(Z = Z, P = P, D = D, Theta = Theta, A = A, n = n, K = K, R = R))
}

###################################################################

gen_data = function(n = 100, K = 1, R = pi / 2, Z = NULL, seed = NULL, seed_Z = NULL){
	# -----------------------------------------------------------------------------
    # Description: Generate data for simulation study with 2-dim hyperbolic model
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # n         - Number of nodes
    # K         - Curvature parameter
    # R         - Radius
    # seed      - Random seed for generating edges
    # seed_Z    - Random seed for generating latent positions

    # -----------------------------------------------------------------------------
    # Output: 
    # Z         - Matrix of latent positions
    # P         - Matrix of linkage probability
    # D         - Matrix of pairwise distance
    # Theta     - Matrix of theta parameter 
    # A         - Adjacency matrix
    # n         - Number of nodes
    # K         - Curvature parameter
    # R         - Radius
    # -----------------------------------------------------------------------------
	if (is.null(Z)){
		set.seed(seed_Z)
		Z = gen_Z(n = n, K = 1, R = R * sqrt(K))
	}
	DPT = gen_P(n = n, K = K, R = R, Z = Z, alpha = alpha)
	P = DPT$P
	D = DPT$D
	Theta = DPT$Theta
	set.seed(seed)
	A = gen_network(n = n, K = K, R = R, P = P, seed = seed)
	return(list(Z = Z, P = P, D = D, Theta = Theta, A = A, n = n, K = K, R = R))
}

###################################################################

gen_network = function(n = 100, K = 1, R = pi / 2, P = NULL, seed = NULL){
	# -----------------------------------------------------------------------------
    # Description: Generate adjacency matrix for simulation study.
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # n         - Number of nodes
    # K         - Curvature parameter
    # R         - Radius
    # P         - Matrix of linkage probability
    # seed      - Random seed for generating edges 

    # -----------------------------------------------------------------------------
    # Output: 
    # A         - Adjacency matrix
    # -----------------------------------------------------------------------------

	if(is.null(P)){
		P = gen_P(n = n, K = K, R = R)$P
	}
	set.seed(seed)
	upper.index <- which(upper.tri(P))
    upper.p <- P[upper.index]
    upper.u <- runif(n = length(upper.p))
    upper.A <- rep(0, length(upper.p))
    upper.A[upper.u < upper.p] <- 1
    A <- matrix(0, n, n)
    A[upper.index] <- upper.A
    A <- A + t(A)
    return(A)
}

###################################################################

gen_P = function(n = 100, K = 1, R = pi / 2, Z = NULL){
	# -----------------------------------------------------------------------------
    # Description: Generate linkage probability matrix for simulation study.
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # n         - Number of nodes
    # K         - Curvature parameter
    # R         - Radius
    # Z         - Matrix of latent positions

    # -----------------------------------------------------------------------------
    # Output: 
    # P         - Matrix of linkage probability
    # D         - Matrix of pairwise distance
    # Theta     - Matrix of theta parameter 
    # -----------------------------------------------------------------------------

	if(is.null(Z)){
		Z = gen_Z(n = n, K = 1, R = sqrt(K) * R)
	}
	if(is.null(alpha)){
		alpha = gen_alpha(n = n)
	}
	d = dim(Z)[2]
	L = diag(c(rep(1, d - 1), -1))
	vec_1 = rep(1, n)
	D = acosh(-Z %*% L %*% t(Z) + diag(rep(1, n))) / sqrt(K)
	diag(D) = rep(0, n)
	Theta = - D

	# generate P with different link function
	# default: logistic(x), then 0 < P < 0.5
	# P = plogis(Theta)
	# modified: logistic(x/2) such that 0 < P < 1
	P = 2 * plogis(Theta)
	# modified: exp(x) such that 0 < P < 1
	# P = exp(Theta)
	diag(P) = rep(0, n)

	return(list(P = P, D = D, Theta = Theta))
}

###################################################################

gen_Z = function(n = 100, K = 1, R = pi / 2){
	# -----------------------------------------------------------------------------
    # Description: Generate latent position uniformly from a hyperbolic disk for simulation study.
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # n         - Number of nodes
    # K         - Curvature parameter
    # R         - Radius

    # -----------------------------------------------------------------------------
    # Output: 
    # (Z1, Z2, Z3)  - Coordinates of latent positions
    # -----------------------------------------------------------------------------

	theta = runif(n, 0, 2 * pi)
	r = runif_h(n, R = R, K = K)
	Z3 = cosh(sqrt(K) * r)
	radius = sqrt(Z3^2 - 1)
	Z1 = radius * cos(theta)
	Z2 = radius * sin(theta)
	return(cbind(Z1 = Z1, Z2 = Z2, Z3 = Z3))
}

###################################################################

gen_data_eu = function(n = 100, R = pi / 2, D = 2, Z = NULL, seed = NULL, seed_Z = 1234){
	# -----------------------------------------------------------------------------
    # Description: Generate data from D-dim Euclidean model for simulation study.
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # n         - Number of nodes
    # R         - Radius
    # D         - Dimension of latent space
    # Z         - Matrix of latent positions
    # seed      - Random seed for generating edges
    # seed_Z    - Random seed for generating latent positions 

    # -----------------------------------------------------------------------------
    # Output: 
    # Z         - Matrix of latent positions
    # P         - Matrix of linkage probability
    # D         - Matrix of pairwise distance
    # Theta     - Matrix of theta parameter 
    # A         - Adjacency matrix
    # n         - Number of nodes
    # K         - Curvature parameter
    # R         - Radius
    # -----------------------------------------------------------------------------

	# generate Z
	if (is.null(Z)){
		set.seed(seed_Z)
		Z = uniformly::runif_in_pball(n, d = D, p = 2, r = R)
	}
	D = as.matrix(dist(Z, method = "euclidean"))
	P = 2 * plogis(-D)
	diag(P) = rep(0, n)

	# generate A
	set.seed(seed)
	A = gen_network(n, P = P)
	return(list(Z = Z, P = P, D = D, Theta = -D, A = A, n = n, K = 0, R = R))

}

###################################################################

gen_tree_network = function(h = 1, m = 4){
	# -----------------------------------------------------------------------------
    # Description: Generate tree-structured network for simulation study.
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # h         - Depth of tree
    # m         - Number of leaves per depth  

    # -----------------------------------------------------------------------------
    # Output: 
    # A         - Adjacency matrix
    # g         - igraph data type for graph

    # -----------------------------------------------------------------------------
	n = (m ** (h + 1) - 1) / (m - 1)
	g = make_tree(n = n, children = m, mode = "undirected")
	A = as_adj(g, type = "both", sparse = FALSE)
	return(list(A = A, g = g))
}










