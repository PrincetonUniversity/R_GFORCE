# project onto constraint set

project_C_perpendicular <- function(Z,S,C,c = 1){
    Z <- (Z + t(Z)) / 2
    S <- (S + t(S)) / 2
    C <- (C + t(C)) / 2
    k1 <- c^2 / (c^2 + 1)
    k2 <- c / (c^2 + 1)
    k3 <- 1 / c
    p <- dim(Z)[1]

    # Solve Mv <- b. v <- [y_a y_T lambda]
    M <- matrix(rep(1,(p+2)^2), ncol = (p+2))
    b <- matrix(rep(0,p+2),ncol = 1)

    #precompute some values
    TC <- sum(diag(C))

    #construct linear system
    M <- M + p*diag(p+2)
    M[1:p,p+2] <- rowSums(C)
    b[1:p] <- rowSums(Z) + k3*rowSums(S)

    M[p+1,1:p] <- 2*rep(1,p)
    M[p+1,p+1] <- p
    M[p+1,p+2] <- TC
    b[p+1] <- sum(diag(Z)) + k3*sum(diag(S))

    M[p+2,1:p] <- 2*colSums(C)
    M[p+2,p+1] <- TC
    M[p+2,p+2] <- sum(C*C)
    b[p+2] <- sum(C*Z) + k3*sum(C*S)

    #compute dual variable values
    v <- qr.solve(M,b)

    #construct projection
    Z_proj <- k1*Z + k2*S - k1*diag(p)*v[p+1] - k1*C*v[p+2]

    j <- matrix(rep(1,p),ncol = 1)
    for(a in 1:p){
        ea <- matrix(rep(0,p),ncol=1)
        ea[a] <- 1
        Ra <- ea%*%t(j) + j%*%t(ea)
        Z_proj <- Z_proj - k1*Ra*v[a]
    }

    res <- NULL
    res$Z_proj <- Z_proj
    res$v <- v
    return(res)
}