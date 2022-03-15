gen_tloc = function(nt, nc){
  function(i){
    (nc*(i-1)+1):(nc*i)
  }
}

gen_cloc = function(nt, nc){
  function(i){
    (0:(nt-1))*nc + i
  }
}

gen_ploc = function(nt, nc){
  tloc = gen_tloc(nt, nc)
  cloc = gen_cloc(nt, nc)
  ploc = function(i, j){
    intersect(tloc(i), cloc(j))
  }
}

approx_int_solution = function(nt, nc, kt, kc, coef, cplex_control=list()){
  tloc = gen_tloc(nt, nc)
  cloc = gen_cloc(nt, nc)

  a2 = matrix(0, nt, nt*nc)
  for (i in 1:nt){
    a2[i, tloc(i)] = 1
  }
  b2 = rep(kt, nt)
  s2 = rep("L", nt)
  
  a3 = matrix(0, nc, nt*nc)
  for (i in 1:nc){
    a3[i, cloc(i)] = 1
  }
  b3 = rep(kc, nc)
  s3 = rep("L", nc)
  
  Amat = rbind(a2, a3)
  bvec = c( b2, b3)
  sense = c(s2, s3)
  cvec = coef
  ub = Inf
  vtype = rep("B", nt*nc)
  
  Rcplex2::Rcplex(cvec, Amat, bvec, sense=sense, ub=ub,
                    objsense='max', vtype = vtype, control=cplex_control)
}

approx_int_solution2 = function(args, coef, cplex_control=list()){
  Rcplex2::Rcplex(Amat=args$Amat,
                        bvec=args$bvec, cvec=coef, ub=args$ub,
                        vtype = "B", n=1, sense=args$sense,
                        objsense='max', control = cplex_control)
}


flip_constraints = function(nt, nc, constraints){
  cloc=gen_cloc(nt, nc)
  flipped_constraints = list()
  for (con in names(constraints)){
    flipped_constraints[[con]] = constraints[[con]]
    if(!is.null(constraints[['con']]$Amat))
      flipped_constraints[[con]] = constraints[[con]]$Amat[, c(t(sapply(1:nc, cloc)))]
    if(!is.null(constraints[[con]]$ub))
      flipped_constraints[[con]]$ub = constraints[[con]]$ub[c(t(sapply(1:nc, cloc)))]
  }
  flipped_constraints
}


cohen = function(tau, Y1, Y0, val){
  n1 = length(Y1)
  n0 = length(Y0)
  mean(Y1 + tau - Y0)/sqrt(((n1 - 1) * var(Y1 + tau) + (n0-1) * var(Y0)) /
                           (n1 + n0 - 2)) - val
}


cohensd  <- function(Y1, Y0){
    n1  <- length(Y1)
    n0  <- length(Y0)
    mu  <- mean(Y1)-mean(Y0)
    s  <- sqrt(((n1-1)*var(Y1) +(n0-1)*var(Y0))/(n1 + n0 - 2))
    mu/s
}

find_cohensd  <- function(beta, W, X, eps1, eps0, dval){
    beta1 <- beta[1:ncol(W)]
    beta0  <- beta[(ncol(W)+1):length(beta)]
    (cohensd(W %*% beta1 + X %*% beta0 + eps1, X%*% beta0 + eps0) - dval)^2
}

gen_from_covariates  <-  function(confounders, moderators,
                                  prop_treated, te_size, var_prop, binary=FALSE){
    W  <-  moderators
    X  <-  confounders
    n  <- nrow(X)

    coefs  <- optim(runif(ncol(W) + ncol(X), -1, 1), find_cohensd, W=W, X=X,
                    eps1=0, eps0=0,
                    dval=te_size, control=list(maxit=10000), method='BFGS')$par

    beta1  <- coefs[1:ncol(W)]
    beta0  <- coefs[(ncol(W)+1):length(coefs)]

    mu0  <- X%*%beta0
    mu1  <- W%*%beta1 + X%*%beta0
    
    D  =  rep(0, n)
    e <- 1/(1+exp(-mu0))
    D[sample(1:n, size=floor(n * prop_treated), replace=F, prob=e)] = 1

    v1 <- abs((W %*% beta1 + X %*% beta0)) * var_prop
    eps1 <- sapply(v1, function(s) {rnorm(1, 0, s)})
    v0 <- abs(X%*% beta0) * var_prop
    eps0 <- sapply(v0, function(s) rnorm(1, 0, s))
    
    Y1 <- mu1 + eps1
    Y0 <- mu0 + eps0

    if(binary){
      Y1 <- Y1 > 0
      Y0 <- Y0 > 0
    }

    Y <- Y1 * D + Y0 * (1-D)
    data.frame(Y1, Y0, Y, D, e)
}

gen_linear  <- function(n, p, prop_treated, te_size,
                        n_moderators, var_prop, binary){
    X  <-  matrix(runif(n*p, -1, 1), n, p)
    W  <-  matrix(1, n, 1)
    if(n_moderators > 0)
        W  <- cbind(W, X[, 1:n_moderators])
    cbind(gen_from_covariates(X, W, prop_treated, te_size, var_prop, binary), X)
}

gen_complex <- function(n, p, prop_treated, te_size,
                        n_moderators, var_prop, binary){
    X  <-  matrix(runif(n*p, -1, 1), n, p)
    W  <-  matrix(1, n, 1)
    if(n_moderators > 0)
        W  <- cbind(W, X[, 1:n_moderators])
    X  <-  cbind(X, X^2)
    X  <-  cbind(X, sin(X))
    X  <-  cbind(X, X[, 1] > 0)
    cbind(gen_from_covariates(X, W, prop_treated, te_size, var_prop, binary), X[, 1:p])
}

gen_stratified <- function(n, p, prop_treated, te_size, n_moderators, var_prop, binary){
  X <- matrix(0, n, p)
  ns <- floor(n/p)
  for (s in 1:p){
    X[(1 + ns*(s-1)):(ns * s), s] <- 1
  }
  W <- matrix(1, n, 1)
  if(n_moderators > 0)
    W  <- cbind(W, X[, 1:n_moderators])
  cbind(gen_from_covariates(X, W, prop_treated, te_size, var_prop, binary), X)
}

McN_exact_p <- function(b, c){
  if (b > c){
    return(1-pbinom(b, b+c, 0.5))
  }else{
    return(pbinom(b, b+c, 0.5))
  }
}