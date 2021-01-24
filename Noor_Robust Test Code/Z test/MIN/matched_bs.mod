## AMPL model for matched pair design
## written by Noor on 22 Sept 2014, edited on 25 sept 2014, copyright Md. Noor-Alam, again edited on jan 02 with revised formulation
## PAL, MIT


set C ={"season", "yr", "workday", "temp", "hum", "windsp"};
# set of all covariates
set n_set={35,45,55,65,75,85,88,90,92,94};
# set of all n

param n default 35;
# total number of pairs
param S;
# sigma temp
param T1=247;
param T2=463;

param NT {i in 1..T1, j in 1..T2} default 0;
param NH {i in 1..T1, j in 1..T2} default 0;
param NW {i in 1..T1, j in 1..T2} default 0;

param DD{i in 1..T1, j in 1..T2} default 0;
param D{i in 1..T1, j in 1..T2} default 0;

param CT {i in 1..T1, k in C};
param CC {j in 1..T2, k in C};

param OT {i in 1..T1}integer;
param OC {j in 1..T2}integer;
param b_l default Infinity;

var a{i in 1..T1,j in 1..T2} integer >=0,<=1 ;


minimize obj1: sum {i in 1..T1,j in 1..T2} ((OT[i]-OC[j])*a[i,j]);

# objective function


subject to constraint1: sum {i in 1..T1,j in 1..T2} ((OT[i]-OC[j])^2)*a[i,j]<=b_l;

subject to max_pair1: sum {i in 1..T1,j in 1..T2}a[i,j]<=n;
subject to max_pair2: sum {i in 1..T1,j in 1..T2}a[i,j]>=n;

subject to avoid_repeat1 {j in 1..T2}: sum {i in 1..T1}a[i,j]<=1;
subject to avoid_repeat2 {i in 1..T1}: sum {j in 1..T2}a[i,j]<=1;

subject to pair_cons {i in 1..T1, j in 1..T2}: a[i,j]<=D[i,j];




