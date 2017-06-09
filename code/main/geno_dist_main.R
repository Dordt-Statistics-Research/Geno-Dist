generate_factor_matrix <- function(nrow, ncol, prob=c(0.25,0.5,0.25)) {
    x <- sample(c(-1,0,1), size=ncol*nrow, replace=T, prob=prob)
    matrix(x,ncol=ncol)
}

generate_parameter_matrix <- function(nrow, ncol, min=-5, max=5) {
    x <- sample(min:max, size=ncol*nrow, replace=T)
    matrix(x,ncol=ncol)
}

generate_default_matrix <- function(nrow, ncol, min=-5, max=5) {
    x <- sample(min:max, size=ncol, replace=T)
    matrix(rep(x,nrow), ncol=ncol, byrow=T)
}

generate_observed_matrix <- function(factor_matrix, parameter_matrix, default_matrix) {
    factor_matrix %*% t(parameter_matrix) + default_matrix
}

partial_derivative <- function(factor_matrix, observation_matrix, pd, i, j) {
    F <- factor_matrix
    Q <- observation_matrix

    F <- cbind(F, 1)
    step1 <- apply(F, 1, function(row){sum(row*pd)})
    step2 <- step1 - Q[,i]
    step3 <- step2 * F[,j]
    2 * sum(step3)
}

gradient <- function(factor_matrix, observation_matrix, pd, i) {
    ## browser()
    F <- factor_matrix
    Q <- observation_matrix

    p <- sapply(seq(ncol(F)), function(k){partial_derivative(F, Q, pd, i, k)})
    d <- partial_derivative(F, Q, pd, i, ncol(F)+1)

    c(p,d)
}

batch_gradient <- function(factor_matrix, observation_matrix, pd, i, vars) {
    
    F <- factor_matrix
    Q <- observation_matrix

    pd <- sapply(vars, function(k){partial_derivative(F, Q, pd, i, k)})

    solution <- rep(0,length(pd))
    solution[vars] <- pd
    solution
}

gradient_descent <- function(factor_matrix, observation_matrix, pd, i, batch=0) {
    min_alpha <- 1e-100
    max_iter <- 10000
    alpha <- 1
    alpha_div <- 2
    ready_to_break <- FALSE
    
    F <- factor_matrix
    Q <- observation_matrix

    g <- function(pd) {
        F1 <- cbind(F, 1)
        step1 <- apply(F1, 1, function(row){sum(pd*row)})
        step2 <- step1 - Q[,i]
        sum(step2 ^ 2)
    }

    if (batch == 0) {
        
    } 
    
    s <- 1
    while(s <= max_iter) {
        ## browser()
        pd_new <- pd - alpha * gradient(F, Q, pd, i)
        if (g(pd_new) < g(pd)) {
            pd <- pd_new
        } else {
            alpha <- alpha / alpha_div   
            if (alpha <= min_alpha) {
                if (ready_to_break) break
                ready_to_break <- T
            }
            s <- s - 1
        }
        s <- s + 1
    }
    ## print(alpha)
    ## print(s)
    ## print(paste("pd_new:",pd_new))
    print(F %*% pd[1:length(pd)-1] + D[1,1] - Q[,1])
    browser()
    pd
}


## F <- generate_factor_matrix(m,n); P <- generate_parameter_matrix(5,n); D <- generate_default_matrix(m,5); Q <- generate_observed_matrix(F,P,D);
