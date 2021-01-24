## AMPL Model for matched pair design with McNemar's test

set C1;
# set of all covariates

param n default 30;
#total number of discordant pairs

param T1;
# total treatment units

param T2;
# total control units

param D{i in 1..T1, j in 1..T2};
#pre-processed binary parameter, which takes a values of 1 if the covariates are close enough to each other, 0 otherwise

param NT {i in 1..T1, j in 1..T2};
param NH {i in 1..T1, j in 1..T2};
param NW {i in 1..T1, j in 1..T2};
param NB {i in 1..T1, j in 1..T2};
param DD{i in 1..T1, j in 1..T2};
# above parameters are required in the run file to calculate D[i,j] from the covariate matrix


param CT {i in 1..T1, k in C1};
#parameters for treatment covariates

param CC {j in 1..T2, k in C1};
#parameters for control covariates

var B;
var C;
#variable for total number of discordant pairs

param OT {i in 1..T1}integer;
#parameters for treatment outcomes

param OC {j in 1..T2}integer;
#parameters for control outcomes

param Q;

var b{i in 1..T1,j in 1..T2} integer >=0,<=1 ;
var c{i in 1..T1,j in 1..T2} integer >=0,<=1 ;

#variables to define discordant pair


var a{i in 1..T1,j in 1..T2} integer >=0,<=1 ;
# if i and j are in same pair then a[i,j]=1, 0 otherwise, also serve for constraint (20)


maximize obj1: (B-C-1)/sqrt(n);
minimize obj2: (B-C-1)/sqrt(n);

# objective function- we can also minimize by replacing maximize as minimize

subject to defb100 {i in 1..T1, j in 1..T2}: a[i,j]*OC[j]*(1-OT[i])>=b[i,j];
subject to defb101 {i in 1..T1, j in 1..T2}: a[i,j]*OC[j]*(1-OT[i])<=b[i,j];

#both of the above constraints represent equation (12) in the paper

subject to defc100 {i in 1..T1, j in 1..T2}: a[i,j]*OT[i]*(1-OC[j])>=c[i,j];
subject to defc101 {i in 1..T1, j in 1..T2}: a[i,j]*OT[i]*(1-OC[j])<=c[i,j];

#both of the above constraints represent equation (13) in the paper

subject to sum1:sum {i in 1..T1,j in 1..T2}b[i,j]>=B;
subject to sum11:sum {i in 1..T1,j in 1..T2}b[i,j]<=B;

#both of the above constraints represent equation (14) in the paper

subject to sum2:sum {i in 1..T1,j in 1..T2}c[i,j]>=C;
subject to sum22:sum {i in 1..T1,j in 1..T2}c[i,j]<=C;

#both of the above constraints represent equation (15) in the paper

subject to max_pair: B+C=n;
#above constraint represents equation (16) in the paper

subject to avoid_repeat2 {i in 1..T1}: sum {j in 1..T2}a[i,j]<=1;
#above constraint represents equation (17) in the paper

subject to avoid_repeat1 {j in 1..T2}: sum {i in 1..T1}a[i,j]<=1;
#above constraint represents equation (18) in the paper


subject to pair_cons {i in 1..T1, j in 1..T2}: a[i,j]<=D[i,j];
#above constraint represents equation (19) in the paper






