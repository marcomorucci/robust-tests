library(openxlsx)
library(designmatch)
library(Matching)
library(xtable)

source("../../robust_test_R/Z_test_Rcplex.R")
source("../../robust_test_R/utils.R")
source("../../constraints_Rcplex.R")
source("../../McN_Rcplex.R")

t_max = 60*5
solver = "cplex"
approximate = 1
solver = list(name = solver, t_max = t_max, approximate = approximate,
              round_cplex = 0, trace_cplex = 0)
cplex_control = list(trace = 0, round = 0, tilim = 60*5, threads=1)
kt = 1
kc = 1

eps = 1e-1
approximate = TRUE
tols = seq(0.01, 1, by=0.01)


#######################################################
## GLOW
#######################################################
glow <- read.xlsx('../data/GLOW.xlsx')
glow$TREATED <- 1-glow$TREATED

X <- glow[, 2:5]
Y <- glow$FRACTURE
D <- glow$TREATED

nt = sum(D)    
nc = sum(1-D)
Xt = X[D==1, ]
Xc = X[D==0, ]
Yt = Y[D==1]
Yc = Y[D==0]
m = nt

Bpairs = as.numeric(sapply(Yt, "*", 1-Yc))
Cpairs = as.numeric(sapply(1-Yt, "*", Yc))

restmat <- NULL
for (i in which(glow$TREATED==1)){
    for(j in which(glow$TREATED==0)){
        if(any(abs(X[i, ] - X[j, ]) > 6)){
            restmat <- rbind(restmat, c(i, j, -1))
        }
    }
}

res <- NULL
# Naive
res = c(res, naive=min(1, 2 * McN_exact_p(sum(Bpairs), sum(Bpairs) + sum(Cpairs))))

## Pscore
psmod = glm(TREATED ~  AGE + WEIGHT + HEIGHT + BMI, 
            data=glow, family=binomial(link='logit'))
ps = Match(Y, D, as.matrix(psmod$fitted.values), caliper=0.1, restrict=restmat,
            estimand = 'ATC', BiasAdjust = F, exact = F, replace=F)
B <- sum(Y[ps$index.treated] * (1-Y[ps$index.control]))
C <- sum(Y[ps$index.control] * (1-Y[ps$index.treated]))
res = c(res, pscore=min(1, 2 * McN_exact_p(B, C)))

## L2
l2 = Match(Y, D, X, estimand = 'ATC', BiasAdjust = F, caliper=1,
             exact = F, replace=F, restrict=restmat)
B <- sum(Y[l2$index.treated] * (1-Y[l2$index.control]))
C <- sum(Y[l2$index.control] * (1-Y[l2$index.treated]))
res = c(res, l2=min(1, 2 * McN_exact_p(B, C)))

## Optimal with moment balance
ord = order(D, decreasing=T)
dmat = distmat(D[ord], X[ord, ])

## Optimal with moment balance
ord = order(D, decreasing=T)
dmat = distmat(D[ord], X[ord, ])

for(i in 1:nrow(dmat)){
    for(j in 1:ncol(dmat)){
        if(any(abs(X[i, ] - X[j, ]) > 6))
            dmat[i, j] <- 1000
    }
}

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
res_glow = c(res, optimal=min(1, 2 * McN_exact_p(B, C)))



#######################################################
## Bike Share
#######################################################

bikeshare <- read.xlsx('../data/BikeShare.xlsx')

X <- bikeshare[, 1:6]
Y <- bikeshare$outcome
D <- bikeshare$treated

nt = sum(D)    
nc = sum(1-D)
Xt = X[D==1, ]
Xc = X[D==0, ]
Yt = Y[D==1]
Yc = Y[D==0]
m = nt

restmat <- NULL
for (i in which(bikeshare$treated==1)){
    for(j in which(bikeshare$treated==0)){
        if(X[i, 'season'] != X[j, 'season'] ||
           X[i, 'yr'] != X[j, 'yr'] ||
           X[i, 'workday'] != X[j, 'workday'] ||
           abs(X[i, 'temp']  - X[j, 'temp']) > 2 ||
           abs(X[i, 'hum'] - X[j, 'hum']) > 6 ||
           abs(X[i, 'windsp'] - X[j, 'windsp'] > 6)){
               restmat <- rbind(restmat, c(i, j, -1))
           }
    }
}

res <- NULL

# Naive
naive_res= t.test(Y[D==1], Y[D==0], alternative = "two.sided")
res = c(res, naive=naive_res$p.value )

## Pscore
psmod <- glm(treated ~ season + yr + workday + temp + hum + windsp, 
                   data=bikeshare,  family=binomial(link='logit'))
ps = Match(Y, D, psmod$fitted.values, restrict=restmat, 
            estimand = 'ATT', M = 1, BiasAdjust = F, exact = F, replace=F)
ps_test <- t.test(Y[ps$index.treated], Y[ps$index.control], alternative='two.sided')

res = c(res, pscore=ps_test$p.value)

## L2
l2 <- Match(Y, D, X, restrict=restmat, 
            estimand = 'ATT', M = 1, BiasAdjust = F, exact = F, replace=F)
l2_test <- t.test(x = Y[l2$index.treated], y=Y[l2$index.control], samedata=FALSE, 
                    alternative="two.sided", bootse=T, bootn=500)
res = c(res, l2=l2_test$p.value)


## Optimal with moment balance
ord = order(D, decreasing=T)

dmat = distmat(D[ord], X[ord, ])
for (i in 1:nrow(dmat)){
    for(j in 1:ncol(dmat)){
        if(X[ord, ][i, 'season'] != X[ord, ][j, 'season'] ||
           X[ord, ][i, 'yr'] != X[ord, ][j, 'yr'] ||
           X[ord, ][i, 'workday'] != X[ord, ][j, 'workday'] ||
           abs(X[ord, ][i, 'temp']  - X[ord, ][j, 'temp']) > 2 ||
           abs(X[ord, ][i, 'hum'] - X[ord, ][j, 'hum']) > 6 ||
           abs(X[ord, ][i, 'windsp'] - X[ord, ][j, 'windsp'] > 6)){
               dmat[i, j] <- 1000
           }
    }
}


# Moment balance: constrain differences in means
# to be at most .05 standard deviations apart
mom_covs = X[ord, ]
mom_tols = round(absstddif(mom_covs, D[ord], .05), 2) 
mom = list(covs = mom_covs, tols = mom_tols)

asd <- capture.output(out <- bmatch(t_ind=D[ord], dist_mat=dmat,
                                    total_groups=nt, subset_weight = NULL,
                                    solver=solver, mom=mom))

opt_test <- t.test(Y[ord][out$t_id], Y[ord][out$c_id], alternative="two.sided")
res_bikeshare = c(res, optimal=t.test(Y[ord][out$t_id], Y[ord][out$c_id],
                            alternative="two.sided")$p.value)

xtable(rbind(res_glow, res_bikeshare))