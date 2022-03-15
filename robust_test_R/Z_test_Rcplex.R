

cplex_Z_problem = function(Yt, Yc, m, kt=1, kc=1, b=Inf, 
                           constraints=list(), approximate=FALSE){
  
  nt = length(Yt)
  nc = length(Yc)
  
  tloc = gen_tloc(nt, nc)
  cloc = gen_cloc(nt, nc)
  ploc = gen_ploc(nt, nc)
  
  cvec = as.numeric(sapply(Yt, "-", Yc))
  
  a1 = matrix(1, 1, nt*nc)
  b1 = m
  s1 = "E"
  
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
  
  a4 = matrix(as.numeric(sapply(Yt, function(y) (y - Yc)**2 )), 1, nt*nc)
  b4 = b
  s4 = "L"
  
  if(approximate)
    vtype="C"
  else
    vtype="B"
  
  args = list(cvec=cvec,
              Amat=rbind(a1, a2, a3, a4),
              bvec=c(b1, b2, b3, b4), 
              sense=c(s1, s2, s3, s4), 
              vtype=vtype,
              ub=1)
  
  for(con in constraints){
    if(!is.null(con$Amat))
      args$Amat = rbind(args$Amat, con$Amat)
    if(!is.null(con$bvec))
      args$bvec = c(args$bvec, con$bvec)
    if(!is.null(con$sense))
      args$sense = c(args$sense, con$sense)
    if(!is.null(con$ub))
      args$ub = con$ub
  }
  args
}

compute_bounded_var = function(d, b, m){
  (d * sqrt(m))/(sqrt(b/m - d^2))
}

maximize_Z_positive = function(nt, nc, m, args, B0, L, eps, 
                               cplex_control=list(trace=0,
                                                  round=1, tilim=60*5), 
                               approximate=FALSE, verbose=TRUE, approxmethod=2){
  
  B = B0
  
  iter = 0
  while(TRUE){
    iter = iter + 1
    if(verbose)
      message(paste0("Iteration ", iter, " grid size = ", length(B)-1))
    
    A = list()
    D = rep(NA, length(B))
    Zub = rep(NA, length(B)-1)
    Zlb = rep(NA, length(B)-1)
    
    if(verbose)
      message("Solving MIPs...")
    
    pb = progress::progress_bar$new(total = length(B)-1)
    for (l in 2:length(B)){
      args$bvec[(1+nt+nc+1)] = B[l]
      al = do.call(Rcplex2::Rcplex, args=c(args, objsense='max',
                                           control=list(cplex_control)))
      pb$tick()
      
      if (al$status == 103)
        next
      
      if (B[l-1]/m > (al$obj/m)^2){
        A[[l]] = al
        D[l] = al$obj/m
        Zub[l] = compute_bounded_var(D[l], B[l], m)
        Zlb[l] = compute_bounded_var(D[l], B[l-1], m)
      }
    }
    zlb = max(Zub, na.rm=T)
    if(verbose)
      message(paste0("Current optimum = ", max(Zlb, na.rm=T),
                     " Solution gap = ", 
                     max(Zlb, na.rm=T) - max(Zub, na.rm=T)))
    if((all(is.na(Zlb)) || all(is.na(Zub))) ||
       (max(Zlb, na.rm = T) - max(Zub, na.rm=T) < eps)){
      break
    }
    
    Bnew = c()
    for(l in 2:length(B)){
      if (!is.na(Zlb[l]) && Zlb[l] > zlb){
        bl = seq(B[l-1], B[l], length.out=3)
        Bnew = c(Bnew, bl)
      }
    }
    B = sort(Bnew)
  }
  
  amax = which.max(Zub)
  if(approximate){
    if(verbose)
      message("Approximating integer solution...")
      if(approxmethod==2)
        aopt = approx_int_solution2(args, A[[amax]]$xopt, cplex_control)
      else
        aopt = approx_int_solution(nt, nc, kt, kc, A[[amax]]$xopt, cplex_control)
  }else
    aopt = A[[amax]]
  dopt = args$cvec %*% aopt$xopt
  Zopt = compute_bounded_var(dopt/m, B[amax], m)
  
  return(list(sol = aopt, Z=Zopt))
}


maximize_Z_negative = function(nt, nc, m, args, B0, L, eps, 
                               cplex_control=list(trace=0, round=1,
                                                  tilim=60*5), 
                               approximate=FALSE, verbose=TRUE, approxmethod=2){
  
  B = B0
  
  iter = 0
  while(TRUE){
    iter = iter + 1
    if(verbose)
      message(paste0("Iteration ", iter, " grid size = ", length(B)-1))
    A = list()
    D = rep(NA, length(B))
    Zub = rep(NA, length(B)-1)
    Zlb = rep(NA, length(B)-1)
    
    if(verbose)
      message("Solving MIPs...")
    pb = progress::progress_bar$new(total = length(B)-1)
    for (l in 2:length(B)){
      args$bvec[(1+nt+nc+1)] = B[l]
      al = do.call(Rcplex2::Rcplex,
                   args=c(args, objsense='min', control=list(cplex_control)))
      
      pb$tick()
      
      if (al$status == 103)
        next
      
      if (B[l-1]/m > (al$obj/m)^2){
        A[[l]] = al
        D[l] = al$obj/m
        Zub[l] = compute_bounded_var(D[l], B[l], m)
        Zlb[l] =  compute_bounded_var(D[l], B[l-1], m)
      }
    }
    
    zlb = min(Zlb, na.rm=T)
    if(verbose)
      message(paste0("Current optimum = ", min(Zub, na.rm=T),
                     " Solution gap = ", min(Zlb, na.rm=T) - min(Zub, na.rm=T)))
    
    if((all(is.na(Zlb)) || all(is.na(Zub))) ||
       (min(Zlb, na.rm=T) - min(Zub, na.rm = T)  < eps))
      break
    
    Bnew = c()
    for(l in 2:length(B)){
      if (!is.na(Zub[l]) && Zub[l] < zlb){
        bl = seq(B[l-1], B[l], length.out=3)
        Bnew = c(Bnew, bl)
      }
    }
    B = sort(Bnew)
    
  }
  amax = which.min(Zub)
  if(approximate){
    if(verbose)
      message("Approximating integer solution...")
    if (approxmethod==2)
      aopt = approx_int_solution2(args, A[[amax]]$xopt, cplex_control)
    else 
      aopt = approx_int_solution(nt, nc, kt, kc, A[[amax]]$xopt, cplex_control)
  }else
    aopt = A[[amax]]
  
  dopt = args$cvec %*% aopt$xopt
  Zopt = compute_bounded_var(dopt/m, B[amax], m)
  
  return(list(sol = aopt, Z=-Zopt))
}

maximize_Z_iterative = function(Y, D, m, kt, kc, L, eps, constraints=list(), 
                                cplex_control=list(trace=0, round=1,
                                                   tilim=60*5),
                                approximate=FALSE, verbose=TRUE, approxmethod=2){
  
  
  Yt = Y[D==1]
  Yc = Y[D==0]
  nt = length(Yt)
  nc = length(Yc)
  
  generated_constraints = list()
  for (el in 1:length(constraints)){
    generated_constraints[[el]] = constraints[[el]](D)
  }
  args = cplex_Z_problem(Yt, Yc, m, kt ,kc, Inf,
                         generated_constraints, approximate)
  a0 = do.call(Rcplex2::Rcplex,
               args=c(args, objsense='max', control=list(cplex_control)))
  if(a0$status == 3){ 
    stop("Matching problem infeasible.")
  }
  
  if(approximate){
    if(approxmethod == 2)
      a0 = approx_int_solution2(args, a0$xopt, cplex_control)
    else
      a0 = approx_int_solution(nt, nc, kt, kc, a0$xopt, cplex_control)

    if(a0$status == 103){
      stop("Matching problem infeasible.")
    }
  }
  
  d0 = args$cvec %*% a0$xopt / m
  b0 = matrix(as.numeric(sapply(Yt,
                                function(y) (y - Yc)**2 )), 1, nt*nc) %*% a0$xopt
  B = suppressWarnings(seq(d0^2*m+1e-3, b0, length.out=L))
  if(verbose)
    message(paste0("Maximal TE = ", d0))
  if (d0 >= 0){
    args = cplex_Z_problem(Yt, Yc, m, kt ,kc,
                           Inf, generated_constraints, approximate)
    return(maximize_Z_positive(nt, nc, m, args, B,
                               L, eps, cplex_control, approximate, verbose, approxmethod))
  }
  if(d0 < 0 ){
    generated_constraints = list()
    for (el in 1:length(constraints)){
      generated_constraints[[el]] = constraints[[el]](1-D)
    }
    args = cplex_Z_problem(Yc, Yt, m, kc, kt, Inf,
                           generated_constraints, approximate)
    return(maximize_Z_negative(nt, nc, m,  args, B, L, eps,
                               cplex_control, approximate, verbose, approxmethod))
  }
}

minimize_Z_iterative = function(Y, D, m, kt, kc, L, eps, constraints=list(), 
                                cplex_control=list(trace=0, round=1, tilim=60*5),
                                approximate=FALSE, verbose=TRUE, approxmethod=2){
  
  sol = maximize_Z_iterative(Y, 1-D, m, kc, kt, L, eps,
                             constraints, cplex_control, 
                             approximate, verbose, approxmethod)
  sol$Z = -sol$Z
  return(sol)
}

## TESTING ###
# n=250
# p = 2
# 
# tau = 5
# alpha = 0
# 
# X = matrix(rnorm(n*p, 0, 1), n, p)
# beta = runif(p, 0, 2)
# 
# e = 0.75/(1+exp(-(X%*%beta/log(p))))
# 
# Y1 = alpha + tau + X%*%beta/log(p) #+ rnorm(n)
# Y0 = alpha + X%*%beta/log(p) #+ rnorm(n)
# D = sapply(e, function(p) rbinom(1, 1, p))
# Y = Y1*D + Y0*(1-D)
# 
# 
# # Put treated first
# X = X[order(D, decreasing=T), ]
# D = D[order(D, decreasing=T)]
# Y = Y[order(D, decreasing=T)]
# 
# 
# Yt = Y[D==1]
# Yc = Y[D==0]
# nt = length(Yt)
# nc = length(Yc)
# 
# Xt = matrix(X[D==1, ], nt, p)
# Xc = matrix(X[D==0, ], nc, p)
# dmat = abs(t(sapply(e, "-", e)))
# 
# m = nt
# kt = 1
# kc = 1
# L = 10
# eps = 1e-2
# approximate = TRUE
# 
# cplex_control = list(trace = 0, round = 1, tilim = 60*5)
# 
# constraints = list(mean=constraint_moment(X, 1, 0.48),
#                    caliper=constraint_caliper(dmat, 0.08))
# max_sol = maximize_Z_iterative(Y, D, m, kt, kc, L, eps, constraints, cplex_control, approximate)
# 
# min_sol = minimize_Z_iterative(Y, D, m, kc, kt, L, eps, constraints, cplex_control, approximate)
# 
# robust_p = min(pnorm(max_sol$Z), 1-pnorm(min_sol$Z))
# 
# 
# 
# # args = cplex_Z_problem(Yt, Yc, m, kt ,kc, Inf, constraints, approximate)
# # a0 = Rcplex::Rcplex(Amat=args$Amat, bvec=args$bvec, cvec=args$cvec, ub=args$ub,
# #                      vtype = args$vtype, n=1, sense=args$sense,
# #                      objsense='max', control = cplex_control)
# # a0 = do.call(Rcplex::Rcplex, args=c(args, objsense='max', control=list(cplex_control)))
# # a1 = Rcplex::Rcplex(Amat=args$Amat, bvec=args$bvec, cvec=a0$xopt, ub=args$ub,
# #                     vtype = "B", n=1,
# #                     objsense='max', control = cplex_control)
# # a0 = approx_int_solution(nt, nc, kt, kc, a0$xopt, cplex_control)
# 
# # Distance matrix
# dist_mat_covs = round(dist(X, diag = TRUE, upper = TRUE), 1)
# dist_mat = as.matrix(dist_mat_covs)
# 
# # Subset matching weight
# subset_weight = 1
# 
# # Moment balance: constrain differences in means to be at most .1 standard deviations apart
# mom_covs = X
# mom_tols = rep(0.1, p)
# mom = list(covs = mom_covs, tols = mom_tols)
# 
# # Solver options
# t_max = 60*5
# solver = "cplex"
# approximate = 1
# solver = list(name = solver, t_max = t_max, approximate = approximate, round_cplex = 1,
#               trace_cplex = 1)
# # Match
# out = bmatch(t_ind=D, dist_mat = NULL, subset_weight = subset_weight, n_controls = 1,
#              mom = mom, solver = solver, total_groups=m)
# 
# # Indices of the treated units and matched controls
# id_1 = out$t_id
# id_2 = out$c_id
# 
# 
# cmat = .constraintmatrix(D, 1, m, mom$covs, mom$tols,  
#                          NULL, NULL, NULL, NULL, NULL, NULL, 
#                          NULL, NULL, NULL, NULL, 
#                          NULL, NULL, NULL, NULL, NULL, 
#                          NULL, NULL, NULL, 0)
# 
# 
# ppars = .problemparameters(D, NULL, 1, 1, m, mom$covs, mom$tols,  
#                            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
#                            NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
#                            approximate)
# constraints = list(con=list(Amat=ppars$Amat[(n+1):(length(ppars$bvec)-1), ], 
#                             bvec=ppars$bvec[(n+1):(length(ppars$bvec)-1)], 
#                             sense=ppars$sense[(n+1):(length(ppars$bvec)-1)]))
# 
# a0 = do.call(Rcplex::Rcplex, args=c(args, objsense='max', n=1,
#                                     control = list(list(trace = 1, 
#                                                         round = 1,
#                                                         tilim = t_max))))
# 
# a1 = Rcplex(cvec=ppars$cvec, Amat=ppars$Amat, bvec=ppars$bvec, ub=ppars$ub, objsense = 'min', n = 1,
#             vtype=ppars$vtype, control = list(trace=1, round=0, tilim=t_max))
