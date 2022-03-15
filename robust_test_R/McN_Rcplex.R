
cplex_max_feasible_pairs = function(Yt, Yc, m, kt, kc, constraints=list(), 
                                    approximate=FALSE){
  nt = length(Yt)
  nc = length(Yc)
  
  tloc = gen_tloc(nt, nc)
  cloc = gen_cloc(nt, nc)

  Bpairs = as.numeric(sapply(Yt, "*", 1-Yc))
  Cpairs = as.numeric(sapply(1-Yt, "*", Yc))
  cvec = Bpairs + Cpairs
  
  a2 = matrix(0, nt, nc*nt)
  for(i in 1:nt){
    a2[i, tloc(i)] = 1
  }
  b2 = rep(kt, nt)
  s2 = rep("L", nt)
  
  a3 = matrix(0, nc, nc*nt)
  for(i in 1:nc){
    a3[i, cloc(i)] = 1
  }
  b3 = rep(kc, nc)
  s3 = rep("L", nc)
  
  a4 = matrix(1, 1, nt*nc)
  b4 = m
  s4 = "E"

  vtype <- ifelse(approximate, "C", "B")

  args = list(cvec=cvec, Amat=rbind(a2, a3, a4), bvec=c(b2,b3, b4),
              sense=c(s2, s3, s4), vtype=vtype, lb=0, ub=1)

  
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

cplex_McN_problem = function(Yt, Yc, m, kt, kc, M, constraints=list(), 
                               approximate=FALSE){
  nt = length(Yt)
  nc = length(Yc)
  
  tloc = gen_tloc(nt, nc)
  cloc = gen_cloc(nt, nc)

  Bpairs = as.numeric(sapply(Yt, "*", 1-Yc))
  Cpairs = as.numeric(sapply(1-Yt, "*", Yc))
  cvec = Bpairs - Cpairs
  
  a1 = Bpairs + Cpairs
  b1 = M
  s1 = "E"
  
  a2 = matrix(0, nt, nc*nt)
  for(i in 1:nt){
    a2[i, tloc(i)] = 1
  }
  b2 = rep(kt, nt)
  s2 = rep("L", nt)
  
  a3 = matrix(0, nc, nc*nt)
  for(i in 1:nc){
    a3[i, cloc(i)] = 1
  }
  b3 = rep(kc, nc)
  s3 = rep("L", nc)
  
  a4 = matrix(1, 1, nt*nc)
  b4 = m
  s4 = "E"

  vtype <- ifelse(approximate, "C", "B")

  args = list(cvec=cvec, Amat=rbind(a1, a2, a3, a4), bvec=c(b1,b2,b3, b4),
              sense=c(s1, s2, s3, s4), vtype=vtype, lb=0, ub=1)

  
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

maximize_McN = function(Y, D, m, kt, kc, minM=NULL, maxM=NULL, L=1, constraints=list(), 
                        approximate=FALSE, approxmethod=2,
                        cplex_control=list(trace=0, round=1, tilim=60*5)){
  Yt = Y[D==1]
  Yc = Y[D==0]
  nt = length(Yt)
  nc = length(Yc)

  generated_constraints <- list()  
  if(!length(constraints)==0){
    for (el in 1:length(constraints)){
      generated_constraints[[el]] = constraints[[el]](D)
    }
  }

  if (is.null(minM))
    minM <- 1

  if (is.null(maxM)){
    args <- cplex_max_feasible_pairs(Yt, Yc, m, kt, kc, 
                                      generated_constraints, approximate)
    a0 = do.call(Rcplex2::Rcplex,
               args=c(args, objsense='max', control=list(cplex_control)))

    if(a0$status == 3){ stop("Matching problem infeasible.")}
    if(approximate){
      if(approxmethod == 2)
        a0 = approx_int_solution2(args, a0$xopt, cplex_control)
      else
        a0 = approx_int_solution(nt, nc, kt, kc, a0$xopt, cplex_control)
    if(a0$status == 3){ stop("Matching problem infeasible.")}
    }
    maxM <- floor(a0$obj)
  }
  
  amax <- NULL
  chimax <- NULL
  argmax <- NULL
  Mmax <- NULL
  for(M in seq(minM, maxM, by=L)){
    args = cplex_McN_problem(Yt, Yc, m, kt ,kc, M,
                          generated_constraints, approximate)

    a0 = do.call(Rcplex2::Rcplex,
                args=c(args, objsense='max', control=list(cplex_control)))

    if(is.null(chimax) || ((a0$obj-1)/sqrt(M) > chimax)){
      amax <- a0
      chimax <- (a0$obj-1)/sqrt(M)
      argmax <- args
      Mmax <- M
    }
  }
  if(is.null(amax)){ stop("Matching problem infeasible.")}

  if(approximate){
      if(approxmethod == 2)
        amax = approx_int_solution2(argmax, amax$xopt, cplex_control)
      else
        amax = approx_int_solution(nt, nc, kt, kc, amax$xopt, cplex_control)

      if(amax$status == 103){
        stop("Matching problem infeasible.")
      }
  } 

  Xopt <- (argmax$cvec %*% amax$xopt - 1)/sqrt(Mmax)
  return(list(sol = amax, chi=Xopt))
}


minimize_McN = function(Y, D, m, kt, kc, minM=NULL, maxM=NULL, L=1, constraints=list(), 
                        approximate=FALSE, approxmethod=2,
                        cplex_control=list(trace=0, round=1, tilim=60*5)){
  Yt = Y[D==1]
  Yc = Y[D==0]
  nt = length(Yt)
  nc = length(Yc)

  generated_constraints <- list()  
  if(!length(constraints)==0){
    for (el in 1:length(constraints)){
      generated_constraints[[el]] = constraints[[el]](D)
    }
  }

  if (is.null(minM))
    minM <- 1

  if (is.null(maxM)){
    args <- cplex_max_feasible_pairs(Yt, Yc, m, kt, kc, 
                                      generated_constraints, approximate)
    a0 = do.call(Rcplex2::Rcplex,
               args=c(args, objsense='max', control=list(cplex_control)))

    if(a0$status == 3){ stop("Matching problem infeasible.")}
    if(approximate){
      if(approxmethod == 2)
        a0 = approx_int_solution2(args, a0$xopt, cplex_control)
      else
        a0 = approx_int_solution(nt, nc, kt, kc, a0$xopt, cplex_control)
    if(a0$status == 3){ stop("Matching problem infeasible.")}
    }
    maxM <- floor(a0$obj)
  }

  amax <- NULL
  chimax <- NULL
  argmax <- NULL
  Mmax <- NULL
  for(M in seq(minM, maxM, by=L)){

    args = cplex_McN_problem(Yt, Yc, m, kt ,kc, M,
                          generated_constraints, approximate)

    a0 = do.call(Rcplex2::Rcplex,
                args=c(args, objsense='min', control=list(cplex_control)))

    if(a0$status == 3){ 
      next
    }

    if(is.null(chimax) || ((a0$obj-1)/sqrt(M) < chimax)){
      amax <- a0
      chimax <- (a0$obj-1)/sqrt(M)
      argmax <- args
      Mmax <- M
    }
  }
  if(is.null(amax)){ stop("Matching problem infeasible.")}

  if(approximate){
      if(approxmethod == 2)
        amax = approx_int_solution2(argmax, amax$xopt, cplex_control)
      else
        amax = approx_int_solution(nt, nc, kt, kc, amax$xopt, cplex_control)

      if(amax$status == 103){
        stop("Matching problem infeasible.")
      }
  } 

  Xopt <- (argmax$cvec %*% amax$xopt - 1)/sqrt(Mmax)
  return(list(sol = amax, chi=Xopt))
}
