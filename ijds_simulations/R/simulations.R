library(Matching, quietly = TRUE)
library(designmatch, quietly = TRUE)

source("../../robust_test_R/Z_test_Rcplex.R")
source("../../robust_test_R/utils.R")
source("../../robust_test_R/constraints_Rcplex.R")
source("../../robust_test_R/McN_Rcplex.R")


z_comparison = function(gen_data, n, p, pt, tau, n_moderators, var_prop, L){
  
  res = c(n=n, p=p, pt=pt, tau=tau, v=var_prop)
  print(res)

  dta <- gen_data(n, p, pt, tau, 0, var_prop, FALSE)
  D <- dta$D
  Y <- dta$Y
  e  <- dta$e
  X  <- as.matrix(dta[, 6:ncol(dta)])
  Y1  <- dta$Y1
  Y0  <- dta$Y0
  
  dmat = abs(t(sapply(e, "-", e)))
  nt = sum(D)    
  nc = sum(1-D)
  Xt = X[D==1, ]
  Xc = X[D==0, ]
  Yt = Y[D==1]
  Yc = Y[D==0]
  m = nt
  M = floor(nt * 0.5)

  for(tol in tols){
    tryCatch({
      mean_tol = sqrt(apply(Xt, 2, var)/2 + apply(Xc, 2, var)/2) * tol
      cal_tol = quantile(dmat[D==1, D==0], tol)
     
      constraints = list(mean=constraint_moment(X, m, 1, mean_tol),
                         caliper=constraint_caliper(dmat, cal_tol))
      
      max_sol = maximize_Z_iterative(Y, D, m, kt, kc, L, eps, constraints, 
                                    cplex_control, approximate, FALSE, approxmethod=1)
      
      min_sol = minimize_Z_iterative(Y, D, m, kc, kt, L, eps, constraints, 
                                    cplex_control, approximate, FALSE, approxmethod=1)

      robust_p  <- min(1,2*min(pnorm(max_sol$Z), 1-pnorm(min_sol$Z)))

      res = c(res, robust=robust_p)

      break
    }, error = function(e) {return(NA)})
  }
  
  # True POs.
  true_res = t.test(Y1, Y0, alternative = "two.sided")
  res = c(res, true=true_res$p.value)
  
  
  # Naive
  naive_res= t.test(Y[D==1], Y[D==0], alternative = "two.sided")
  res = c(res, naive=naive_res$p.value )
  
  
  ## Pscore
  ps = glm(D ~ X, family=binomial)
  ps = Match(Y, D, ps$fitted.values, 
             estimand = 'ATT', M = 1, BiasAdjust = F, exact = F, replace=F)
  
  res = c(res, pscore=2 * (1- pnorm(abs(ps$est/ps$se.standard))))
  
  ## L2
  l2 = Match(Y, D, X, 
             estimand = 'ATT', M = 1, BiasAdjust = F, exact = F, replace=F)
  
  res = c(res, l2=min(1, 2 * (1- pnorm(abs(l2$est/l2$se.standard)))))
  
  ## Optimal with moment balance
  ord = order(D, decreasing=T)
  dmat = distmat(D[ord], X[ord, ])
  
  # Moment balance: constrain differences in means
  # to be at most .05 standard deviations apart
  mom_covs = X[ord, ]
  mom_tols = round(absstddif(mom_covs, D[ord], .05), 2)
  mom = list(covs = mom_covs, tols = mom_tols)
  
  asd <- capture.output(out <- bmatch(t_ind=D[ord], dist_mat=dmat,
                                      total_groups=nt, subset_weight = NULL,
                                      solver=solver, mom=mom))
  res = c(res, optimal=t.test(Y[ord][out$t_id], Y[ord][out$c_id],
                              alternative="two.sided")$p.value)
  
  res
}


approx_comparison = function(gen_data, n, p, pt, tau, n_moderators, var_prop, L){
  
  res = c(n=n, p=p, pt=pt, tau=tau, v=var_prop)
  print(res)

  dta <- gen_data(n, p, pt, tau, 0, var_prop, FALSE)
  D <- dta$D
  Y <- dta$Y
  e  <- dta$e
  X  <- as.matrix(dta[, 6:ncol(dta)])
  Y1  <- dta$Y1
  Y0  <- dta$Y0
  
  dmat = abs(t(sapply(e, "-", e)))
  nt = sum(D)    
  nc = sum(1-D)
  Xt = X[D==1, ]
  Xc = X[D==0, ]
  Yt = Y[D==1]
  Yc = Y[D==0]
  m = nt

  for(tol in tols){
    tryCatch({
      #mean_tol = sqrt(apply(Xt, 2, var)/2 + apply(Xc, 2, var)/2) * tol
      cal_tol = quantile(dmat[D==1, D==0], tol)
     
      constraints = list(#mean=constraint_moment(X, m, 1, mean_tol),
                         caliper=constraint_caliper(dmat, cal_tol))
      t0 <- proc.time()[3]
      max_sol = maximize_Z_iterative(Y, D, m, kt, kc, L, eps, constraints, 
                                    cplex_control, FALSE, FALSE)
      
      min_sol = minimize_Z_iterative(Y, D, m, kc, kt, L, eps, constraints, 
                                    cplex_control, FALSE, FALSE)
      mip_time <- proc.time()[3] - t0
      mip_p  <- min(1, 2*min(pnorm(max_sol$Z), 1-pnorm(min_sol$Z)))

      t0 <- proc.time()[3]
      max_sol = maximize_Z_iterative(Y, D, m, kt, kc, L, eps, constraints, 
                                    cplex_control, TRUE, FALSE, 1)
      
      min_sol = minimize_Z_iterative(Y, D, m, kc, kt, L, eps, constraints, 
                                    cplex_control, TRUE, FALSE, 1)
      a1_time <- proc.time()[3] - t0
      a1_p  <- min(1, 2*min(pnorm(max_sol$Z), 1-pnorm(min_sol$Z)))

      t0 <- proc.time()[3]
      max_sol = maximize_Z_iterative(Y, D, m, kt, kc, L, eps, constraints, 
                                    cplex_control, TRUE, FALSE, 2)
      
      min_sol = minimize_Z_iterative(Y, D, m, kc, kt, L, eps, constraints, 
                                    cplex_control, TRUE, FALSE, 2)
      a2_time <- proc.time()[3] - t0
      a2_p  <- min(1, 2*min(pnorm(max_sol$Z), 1-pnorm(min_sol$Z)))


      res = c(res, mip=mip_p, approximation_1=a1_p, approximation_2=a2_p, 
              mip_time=mip_time, a1_time=a1_time, a2_time=a2_time)
      break
    }, error = function(e) {return(NA)})
  }
  res
}


mcn_comparison = function(gen_data, n, p, pt, tau, n_moderators=0, var_prop=NULL, L){
  
  res = c(n=n, p=p, pt=pt, tau=tau, v=var_prop)
  print(res)

  dta <- gen_data(n, p, pt, tau, n_moderators, var_prop, TRUE)
  D <- dta$D
  Y <- dta$Y
  e  <- dta$e
  X  <- as.matrix(dta[, 6:ncol(dta)])
  Y1  <- dta$Y1
  Y0  <- dta$Y0

  dmat = abs(t(sapply(e, "-", e)))
  nt = sum(D)    
  nc = sum(1-D)
  Xt = X[D==1, ]
  Xc = X[D==0, ]
  Yt = Y[D==1]
  Yc = Y[D==0]
  m = nt

  Bpairs = as.numeric(sapply(Yt, "*", 1-Yc))
  Cpairs = as.numeric(sapply(1-Yt, "*", Yc))

  for(tol in tols){
    tryCatch({
      t0 <- proc.time()[3]
      mean_tol = sqrt(apply(Xt, 2, var)/2 + apply(Xc, 2, var)/2) * tol
      cal_tol = quantile(dmat[D==1, D==0], tol)
     
      constraints = list(mean=constraint_moment(X, m, 1, mean_tol),
                         caliper=constraint_caliper(dmat, cal_tol))
      
      max_sol = maximize_McN(Y, D, m, kt, kc, 1, NULL, L, constraints, 
                             approximate, 1, cplex_control)
      
      min_sol = minimize_McN(Y, D, m, kt, kc, 1, NULL, L, constraints, 
                             approximate, 1, cplex_control)
      
      pmax <- pbinom(Bpairs %*% max_sol$sol$xopt, 
                     Bpairs %*% max_sol$sol$xopt + Cpairs %*% max_sol$sol$xopt, 0.5)
      pmin <- 1- pbinom(Bpairs %*% min_sol$sol$xopt, 
                        Bpairs %*% min_sol$sol$xopt + Cpairs %*% min_sol$sol$xopt, 0.5)

      robust_p  <- min(1, 2*min(pmax, pmin))

      res = c(res, robust=robust_p, proc.time()[3] - t0)
     
      break
    }, error = function(e) {return(NA)})
  }
  
  # True POs.
  B <- sum(Y1 * (1-Y0))
  C <- sum(Y0 * (1-Y1))
  res = c(res, true=2*McN_exact_p(B, C))
  
  # Naive
  res = c(res, naive=min(1, 2 * McN_exact_p(sum(Bpairs), sum(Bpairs) + sum(Cpairs))))
  
  ## Pscore
  ps = glm(D ~ X, family=binomial)
  ps = Match(Y, D, ps$fitted.values, 
             estimand = 'ATT', M = 1, BiasAdjust = F, exact = F, replace=F)
  B <- sum(Y[ps$index.treated] * (1-Y[ps$index.control]))
  C <- sum(Y[ps$index.control] * (1-Y[ps$index.treated]))
  res = c(res, pscore=min(1, 2 * McN_exact_p(B, C)))
  
  ## L2
  l2 = Match(Y, D, X, 
             estimand = 'ATT', M = 1, BiasAdjust = F, exact = F, replace=F)
  B <- sum(Y[l2$index.treated] * (1-Y[l2$index.control]))
  C <- sum(Y[l2$index.control] * (1-Y[l2$index.treated]))
  res = c(res, l2=min(1, 2 * McN_exact_p(B, C)))
  
  ## Optimal with moment balance
  ord = order(D, decreasing=T)
  dmat = distmat(D[ord], X[ord, ])
  
  # Moment balance: constrain differences in means
  # to be at most .05 standard deviations apart
  mom_covs = X[ord, ]
  mom_tols = round(absstddif(mom_covs, D[ord], .05), 2)
  mom = list(covs = mom_covs, tols = mom_tols)
  
  asd <- capture.output(out <- bmatch(t_ind=D[ord], dist_mat=dmat,
                                      total_groups=nt, subset_weight = NULL,
                                      solver=solver, mom=mom))

  B <- sum(Y[ord][out$t_id] * (1 - Y[ord][out$c_id]))
  C <- sum(Y[ord][out$c_id] * (1 - Y[ord][out$t_id]))
  res = c(res, optimal=min(1, 2 * McN_exact_p(B, C)))
  
  res
}


mb_comparison = function(gen_data, n, p, pt, tau, n_moderators, var_prop, L){
  
    res = c(n=n, p=p, pt=pt, tau=tau, v=var_prop)
    print(res)

    dta <- gen_data(n, p, pt, tau, 0, var_prop, FALSE)

    D <- dta$D
    Y <- dta$Y
    e  <- dta$e
    X  <- as.matrix(dta[, 6:ncol(dta)])
    Y1  <- dta$Y1
    Y0  <- dta$Y0
    
    dmat = abs(t(sapply(e, "-", e)))
    nt = sum(D)    
    nc = sum(1-D)
    Xt = X[D==1, ]
    Xc = X[D==0, ]
    Yt = Y[D==1]
    Yc = Y[D==0]
    m = sum(D)
 

    for(tol in tols){
        tryCatch({
            mean_tol = sqrt(apply(Xt, 2, var)/2 + apply(Xc, 2, var)/2) * tol
            cal_tol = quantile(dmat[D==1, D==0], tol)
            
            constraints = list(mean=constraint_moment(X, m, 1, mean_tol),
                            caliper=constraint_caliper(dmat, cal_tol))
            
            max_sol = maximize_Z_iterative(Y, D, m, kt, kc, L, eps, constraints, 
                                            cplex_control, approximate, FALSE,
                                            approxmethod=1)
            max_m = which(matrix(max_sol$sol$xopt, nt, nc, byrow=T)==1, arr.ind=T)
            max_te_rt = mean(Yt[max_m[, 1]] - Yc[max_m[, 2]])
            
            min_sol = minimize_Z_iterative(Y, D, m, kc, kt, L, eps, constraints, 
                                            cplex_control, approximate, FALSE,
                                            approxmethod=1)
            min_m = which(matrix(min_sol$sol$xopt, nt, nc, byrow=T)==1, arr.ind=T)
            min_te_rt = mean(Yt[min_m[, 1]] - Yc[min_m[, 2]])
            
            res = c(res, rt_pval=2*min(pnorm(max_sol$Z), 1-pnorm(min_sol$Z)),    
                    rt_max_att=max_te_rt, rt_min_att = min_te_rt)
            
            generated_constraints = list()
            for (el in 1:length(constraints)){
                generated_constraints[[el]] = constraints[[el]](D)
            }
            args = cplex_Z_problem(Yt, Yc, m, kt ,kc, Inf, 
                generated_constraints, approximate)
            max_mb = do.call(Rcplex2::Rcplex,
                args=c(args, objsense='max', control=list(cplex_control)))
            max_mb = approx_int_solution(nt, nc, kt, kc, max_mb$xopt, cplex_control)
            max_m = which(matrix(max_mb$xopt, nt, nc, byrow=T)==1, arr.ind=T)
            max_te_mb = mean(Yt[max_m[, 1]] - Yc[max_m[, 2]])
            
            min_mb = do.call(Rcplex2::Rcplex, 
                args=c(args, objsense='min', control=list(cplex_control)))
            min_mb = approx_int_solution(nt, nc, kt, kc, min_mb$xopt, cplex_control)
            min_m = which(matrix(min_mb$xopt, nt, nc, byrow=T)==1, arr.ind=T)
            min_te_mb = mean(Yt[min_m[, 1]] - Yc[min_m[, 2]])
            
            max_mb_p = t.test(args$cvec[max_mb$xopt >0], alternative = "less")$p.value
            min_mb_p = t.test(args$cvec[min_mb$xopt >0], alternative = "greater")$p.value
            mb_p = min(max_mb_p, min_mb_p)

            res = c(res, mb_pval=mb_p, mb_max_att=max_te_mb, mb_min_att=min_te_mb)
            break
        }, error = function(e) {return(NA)})
    }
  res
}

constraint_comparison = function(gen_data, n, p, pt, tau, n_moderators, var_prop, L){
  
  res = c(n=n, p=p, pt=pt, tau=tau, v=var_prop)
  print(res)

  dta <- gen_data(n, p, pt, tau, 0, var_prop, FALSE)
  D <- dta$D
  Y <- dta$Y
  e  <- dta$e
  X  <- as.matrix(dta[, 6:ncol(dta)])
  Y1  <- dta$Y1
  Y0  <- dta$Y0
  
  dmat = abs(t(sapply(e, "-", e)))
  nt = sum(D)    
  nc = sum(1-D)
  Xt = X[D==1, ]
  Xc = X[D==0, ]
  Yt = Y[D==1]
  Yc = Y[D==0]
  m = sum(D)

  for(tol in tols){
      tryCatch({
        mean_tol = sqrt(apply(Xt, 2, var)/2 + apply(Xc, 2, var)/2) * tol
        var_tol = (apply(Xt, 2, var)/2 + apply(Xc, 2, var)/2) * tol
        mom_constraints = list(mean=constraint_moment(X, m, 1, mean_tol), 
                               variance=constraint_moment(X, m, 2, var_tol))
        t0 = proc.time()[3]
        mom_max_sol = maximize_Z_iterative(Y, D, m, kt, kc, L, eps, mom_constraints, 
                                           cplex_control, approximate, FALSE, 1)
        mom_min_sol = minimize_Z_iterative(Y, D, m, kc, kt, L, eps, mom_constraints, 
                                           cplex_control, approximate, FALSE, 1)
        res = c(res, mom_time=proc.time()[3] - t0)
        res = c(res, mom_pval=2*min(pnorm(mom_max_sol$Z), 1-pnorm(mom_min_sol$Z)))
        break
      }, error=function(e) return(NA))
    }
    
    for(tol in tols){
      tryCatch({
        cal_constraints = list(caliper=constraint_caliper(dmat, 
                                caliper=quantile(dmat[D==1, D==0], tol)))
        t0 = proc.time()[3]
        cal_max_sol = maximize_Z_iterative(Y, D, m, kt, kc, L, eps, cal_constraints, 
                                           cplex_control, approximate, FALSE, 1)
        cal_min_sol = minimize_Z_iterative(Y, D, m, kc, kt, L, eps, cal_constraints, 
                                           cplex_control, approximate, FALSE, 1)
        res = c(res, cal_time=proc.time()[3] - t0)
        res = c(res, cal_pval=2*min(pnorm(cal_max_sol$Z), 1-pnorm(cal_min_sol$Z)))
        break
          }, error=function(e) return(NA))
    }
    
    for(tol in tols){
      tryCatch({
        quantiles = sapply(1:p, function(j) seq(min(X[, j]), max(X[, j]), length.out=5), 
                           simplify = F)
        qua_constraints = list(quantile=constraint_quantile(X, m, quantiles, tol))
        t0 = proc.time()[3]
        qua_max_sol = maximize_Z_iterative(Y, D, m, kt, kc, L, eps, qua_constraints, 
                                           cplex_control, approximate, FALSE, 1)
        qua_min_sol = minimize_Z_iterative(Y, D, m, kc, kt, L, eps, qua_constraints, 
                                           cplex_control, approximate, FALSE, 1)
        res = c(res, qua_time=proc.time()[3] - t0)
        res = c(res, qua_pval=min(pnorm(qua_max_sol$Z), 1-pnorm(qua_min_sol$Z)))
        break
          }, error=function(e) return(NA))
  }
  res
}

run_sim_slurm  <-  function(sim_fcn, data_gen, N, ATE, P, PT, VAR){
    all_res  <- NULL
    for (n in N){
      for(tau in ATE){
        for(p in P){
          for (vr in VAR){
            for(pt in PT){
              L <- floor(n * 0.2)
              sres  <- sim_fcn(data_gen, n, p, pt, tau, 0, vr, L)
              all_res = rbind(all_res, sres)
            }
          }
        }
      }
    }
    all_res
}
