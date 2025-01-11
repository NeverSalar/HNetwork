
runif_h_3d = function(n = 100, R = pi / 2, K = 1) {
    # -----------------------------------------------------------------------------
    # Description: Uniform sampling on 3-dim Poincare ball.
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # n         - Number of nodes
    # R         - Radius
    # K         - Curvature parameter
    # -----------------------------------------------------------------------------
    # Output: 
    # (x, y, z) - Coordinates of latent positions
    # -----------------------------------------------------------------------------

    # Step 1: Sample radius
    u <- runif(n, 0, 1)
    r <- acosh(1 + u * (cosh(sqrt(K) * R) - 1)) / sqrt(K)

    # Step 2: Sample angular coordinates
    v <- runif(n, 0, 1)
    theta <- acos(1 - 2 * v)
    phi <- runif(n, 0, 2 * pi)

    # Step 3: Convert to Cartesian coordinates
    x <- r * sin(theta) * cos(phi)
    y <- r * sin(theta) * sin(phi)
    z <- r * cos(theta)
    return(cbind(x, y, z))
}

###################################################################

runif_h = function(n = 100, R = pi / 2, K = 1){
    # -----------------------------------------------------------------------------
    # Description: Uniform sampling on 2-dim Poincare disk.
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # n         - Number of nodes
    # R         - Sampling radius
    # K         - Curvature parameter
    # -----------------------------------------------------------------------------
    # Output: 
    # (x, y)    - Coordinates of latent positions
    # -----------------------------------------------------------------------------

    # use inverse transformation sampling
    u = runif(n, 0, 1)
    res = acosh(1 + (cosh(sqrt(K) * R) - 1) * u) / sqrt(K)
    return(res)
}

###################################################################

hyperboloid_to_Poincare = function(Z){
    # -----------------------------------------------------------------------------
    # Description: Mapping hyperboloid coordinates to Poincare coordinates.
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # Z         - Matrix of coordinates under hyperboloid model
    # -----------------------------------------------------------------------------
    # Output: 
    # Z_P       - Matrix of coordinates under Poincare model
    # -----------------------------------------------------------------------------
    n = dim(Z)[1]
    d = dim(Z)[2]
    Z_P = c()
    for (i in 1 : n){
        temp_coor = c(Z[i, 1 : (d - 1)]) / (Z[i, d] + 1)
        # temp_coor = c(Z[i, 1], Z[i, 2]) / (Z[i, 3] + 1)
        Z_P = rbind(Z_P, temp_coor)
    }
    return(Z_P)
}

###################################################################

plot_Poincare = function(Z_P, label = NULL, degrees = NULL, main = NULL){
    # -----------------------------------------------------------------------------
    # Description: Visualizing hyperbolic embeddings, supporting coloring nodes with labels and sizing with degrees.
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # Z_P         - Matrix of coordinates under Poincare model
    # label       - Group labels used for coloring nodes
    # degrees     - Node degree
    # main        - Title of plot
    # -----------------------------------------------------------------------------
    # Output: 
    # Plot        - Visualization of latent positions
    # -----------------------------------------------------------------------------
    
    n = dim(Z_P)[1]

    # plot the unit circle
    tt <- seq(0,2 * pi,length.out = 100) # seqence of hundred points from 0 to 2-phi .
    xx <- cos(tt) # x cordinate of circle with the sequence above.
    yy <- sin(tt) # y cordinate of circle with the sequence above.
    
    if (is.null(main)){
        main = "Visualization of Latent Positions with Poincare Model"
    }
    plot(xx, yy, type = "l", xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), xlab = '', ylab = '', main = main)
    par(new = TRUE) # "hold" the plot above.
    # plot the data points
    plot(Z_P[ ,1], Z_P[ ,2], xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), xlab = '', ylab = '', col = label, cex = degrees, lwd = 2) # plot the embeddings on top of the circle.
}

###################################################################

dist_to_Poincare_center = function(x){
    # -----------------------------------------------------------------------------
    # Description: Calculate distance between a point on Poincare disk and disk center.
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # x           - Coordicates of the point
    # -----------------------------------------------------------------------------
    # Output: 
    # d           - Distance to the center
    # -----------------------------------------------------------------------------
    temp = 2 * x ^ 2 / (1 - x ^ 2)
    d = acosh(1 + temp)
    return(d)
}

###################################################################

dist_two_points_on_Poincare_disk = function(a, b){
    # -----------------------------------------------------------------------------
    # Description: Calculate distance between two points on Poincare disk.
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # a, b        - Coordinates of two points
    # -----------------------------------------------------------------------------
    # Output: 
    # d           - Distance between two points
    # -----------------------------------------------------------------------------

    # point a and b, two-dim vectors
    temp = 2 * sum((a - b)^2) / (1 - sum(a^2)) / (1 - sum(b^2))
    d = acosh(1 + temp)
    return(d)
}

###################################################################

iso_trans_x_axis_loid = function(theta, Z){
    # -----------------------------------------------------------------------------
    # Description: Rotate points on hyperboloid by theta with x axis.
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # theta       - Angle of rotation
    # Z           - Matrix of latent positions
    # -----------------------------------------------------------------------------
    # Output: 
    # Z % Ax      - Rotated latent positions
    # -----------------------------------------------------------------------------
    Ax = matrix(c(1, 0, 0, 0, cosh(theta), sinh(theta), 0, sinh(theta), cosh(theta)), ncol = 3)
    return(Z %*% Ax)
}

###################################################################

iso_trans_y_axis_loid = function(theta, Z){
    # -----------------------------------------------------------------------------
    # Description: Rotate points on hyperboloid by theta with y axis.
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # theta       - Angle of rotation
    # Z           - Matrix of latent positions
    # -----------------------------------------------------------------------------
    # Output: 
    # Z % Ay      - Rotated latent positions
    # -----------------------------------------------------------------------------
    Ay = matrix(c(cosh(theta), 0, sinh(theta), 0, 1, 0, sinh(theta), 0, cosh(theta)), ncol = 3)
    return(Z %*% Ay)
}

###################################################################

iso_rot_Poincare = function(theta, Z){
    # -----------------------------------------------------------------------------
    # Description: Rotate points on hyperboloid by theta with z axis, around (0, 0).
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # theta       - Angle of rotation
    # Z           - Matrix of latent positions
    # -----------------------------------------------------------------------------
    # Output: 
    # Z % Az      - Rotated latent positions
    # -----------------------------------------------------------------------------
    Az = matrix(c(cos(theta), sin(theta), 0, - sin(theta), cos(theta), 0, 0, 0, 1), ncol = 3)
    return(Z %*% Az)
}

###################################################################

project_E_to_H = function(Z, d = 3, Rmax = cosh(sqrt(K) * 4)){
    # -----------------------------------------------------------------------------
    # Description: Project Eulidean points on Poincare disk.
    # Date: Jan, 2025
    # -----------------------------------------------------------------------------
    # Input:
    # Z           - Matrix of latent positions
    # d           - Dimension of coordinates (= dimension of hyperbolic space - 1)
    # Rmax        - The maximum distance of projected point to disk center
    # -----------------------------------------------------------------------------
    # Output: 
    # t(res)      - Matrix of projected positions
    # -----------------------------------------------------------------------------
    n = dim(Z)[1]
    Lambda = diag(c(rep(1, d - 1), -1))
    I = diag(rep(1, d))
    res = sapply(1 : n, function(i){
        fun = function(l){
            Zz = (I + l * Lambda) %*% Z[i, ]
            return(t(Zz) %*% Lambda %*% Zz + 1)}
        lambda = nleqslv(c(0), fun, control = list(maxit = 10000))$x
        new_Zi = (I + lambda * Lambda) %*% Z[i, ]
        # restrict all embeddings within a disk with radius 4
        if(new_Zi[d] > Rmax){
            new_Zi = new_Zi / new_Zi[d] * Rmax
            new_Zi[d] = sqrt(1 + sum(new_Zi[1:(d-1)] * new_Zi[1:(d-1)]))
        }
        return(new_Zi)})
    return(t(res))
}

