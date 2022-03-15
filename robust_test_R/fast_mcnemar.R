prune_excess = function(Yt, Yc, maximize=TRUE){
  nt = length(Yt)
  nc = length(Yc)
  if(maximize){
    if(nc > nt)
      return(cbind(treated=Yt, control=sort(Yc)[1:nt]))
    else if (nt == nc)
      return(cbind(treated=Yt, control=Yc))
    else
      return(cbind(treated=sort(Yt, decreasing=T)[1:nc], control=Yc))
  }else{
    pr_vals = prune_excess(Yc, Yt, TRUE)
    return(cbind(treated=pr_vals[, 2], control=pr_vals[, 1])) 
  }
}


compute_max_mcnemar = function(Yts, Ycs){
  N = 0
  TE = 0.0
  SD = 0.0
  for (l in 1:length(Yts)){
    if(length(Yts[[l]])==0 || length(Ycs[[l]])==0)
      next
    pouts = prune_excess(Yts[[l]], Ycs[[l]], TRUE)
    N = N + nrow(pouts)
    TE = TE + sum(pouts[,1]) - sum(pouts[, 2])
  }
  
  if (TE > 0){
    for (l in 1:length(Yts)){
      if(length(Yts[[l]])==0 || length(Ycs[[l]])==0)
        next
      pouts = prune_excess(Yts[[l]], Ycs[[l]], TRUE)
      SD = SD + abs(sum(pouts[,1]) - sum(pouts[, 2]))
    }
  }else{
    for (l in 1:length(Yts)){
      if(length(Yts[[l]])==0 || length(Ycs[[l]])==0)
        next
      pouts = prune_excess(Yts[[l]], Ycs[[l]], TRUE)
      SD = SD + nrow(pouts) - abs(sum(pouts[, 1]) + sum(pouts[, 2]) - nrow(pouts))
    }
  }
  
  return ((TE - 1)/ sqrt(SD + 1))
}



compute_min_mcnemar = function(Yts, Ycs){
  N = 0
  TE = 0.0
  SD = 0.0
  for (l in 1:length(Yts)){
    if(length(Yts[[l]])==0 || length(Ycs[[l]])==0)
      next
    pouts = prune_excess(Yts[[l]], Ycs[[l]], FALSE)
    N = N + nrow(pouts)
    TE = TE + sum(pouts[,1]) - sum(pouts[, 2])
  }
  
  if (TE <= 0){
    for (l in 1:length(Yts)){
      if(length(Yts[[l]])==0 || length(Ycs[[l]])==0)
        next
      pouts = prune_excess(Yts[[l]], Ycs[[l]], FALSE)
      SD = SD + abs(sum(pouts[,1]) - sum(pouts[, 2]))
    }
  }else{
    for (l in 1:length(Yts)){
      if(length(Yts[[l]])==0 || length(Ycs[[l]])==0)
        next
      pouts = prune_excess(Yts[[l]], Ycs[[l]], FALSE)
      SD = SD + nrow(pouts) - abs(sum(pouts[, 1]) + sum(pouts[, 2]) - nrow(pouts))
    }
  }
  
  return ((TE - 1)/ sqrt(SD + 1))
}

