constraint_caliper = function(dist_matrix, caliper){
  add_constraint_caliper = function(D){
    dmat = dist_matrix[D==1, D==0]
    ub = as.numeric(as.numeric(t(dmat)) <= caliper)
    list(ub=ub)
  }
}

constraint_moment = function(X, m, mom, eps){
  add_constraint_moment = function(D){
    Xt = X[D==1, ]
    Xc = X[D==0, ]

    nt = nrow(Xt)
    nc = nrow(Xc)
    p = ncol(Xt)

    if(length(eps) == 1)
      eps = rep(eps, p)

    a = matrix(0, p, nt*nc)
    for (j in 1:p){
      a[j, ] = as.numeric(sapply(Xt[, j]**mom, "-", Xc[, j]**mom))/m
    }
    b = eps
    s = rep("L", p)

    Amat = rbind(a, -a)
    bvec = c(b, b)
    sense = c(s, s)
    list(Amat=Amat, bvec=bvec, sense=sense)
  }
}

constraint_quantile = function(X, m, quantiles, eps){
  
  add_constraint_quantile = function(D){
    Xt = X[D==1, ]
    Xc = X[D==0, ]
    nt = nrow(Xt)
    nc = nrow(Xc)
    p = ncol(Xt)
    
    a = matrix(0, sum(sapply(quantiles, length)), nt*nc)
    i = 1
    for (j in 1:p){
      for(q in quantiles[[j]]){
        a[i, ] = as.numeric(sapply((Xt[, j] <= q), "-", (Xc[, j] <= q)))/m
        i = i + 1
      }
    }
    if (length(eps)==1)
      b = rep(eps, sum(sapply(quantiles, length)))
    else
      b = unlist(eps)
    s = rep("L", sum(sapply(quantiles, length)))
    

    Amat = rbind(a, -a)
    bvec = c(b, b)
    sense = c(s, s)
    list(Amat=Amat, bvec=bvec, sense=sense)
  }
}
