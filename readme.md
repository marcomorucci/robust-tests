# Replication Code for the Paper: A robust approach to quantifying uncertainty in matching problems of causal inference

Link to the paper: https://arxiv.org/abs/1812.02227

Ready to use code for all the methods described in the paper is available in the robust_test_R subfolder of this repository. 

*R platform used*
- R version 4.0.4 (2021-02-15)
- Platform: x86_64-centos-linux-gnu (64-bit)
- Running under: Red Hat Enterprise Linux 8.4 (Ootpa)

*R packages required to run the replication code (version used):* 
- reshape2 (1.4.4), ggplot2 (3.3.5), dplyr (1.0.8), Matching (4.9-11), designmatch (0.3.1), openxlsx (4.2.5), xtable (1.8-4), slam (0.1-50), lattice (0.20-41), MASS (7.3-56), Rcplex (0.3-3), Rgplk (0.6-4)

*Python platform used*
- Python 3.8.4 (v3.8.4:dfa645a65e, Jul 13 2020, 10:45:06), [Clang 6.0 (clang-600.0.57)] on darwin

*Python packages required to run the code (version used):*
- numpy (1.20.3), scipy (1.6.3), pandas (1.2.4) (for Figures 2 and 10)
- numpy (1.20.3), scipy (1.6.3), matplotlib (3.4.2), tabulate (0.8.7) (for Figure 9)

*Software required: AMPL (20200501), Excel, GLPK (4.65), Cplex (20.1.0 - Simulations), Cplex (12.10 - Case studies) *

**Before trying to run any of the files in your own environment please make sure that all the paths point to the correct folders.**

### To reproduce simulation results
    run: ijds_simulations/<fname>.sh to submit the simulation job to a SLURM cluster
    run: make_outputs.R to create plots from the simulation results

### To reproduce results in Figure 9 (appendix) 
    run: distribution_comparison.py

### To reproduce results in Table 4 (appendix)
    run: ijds_simulations/case_studies_comparison.R

### To reproduce results in Case Study 1:
    Load the following files in AMPL: covariates.dat, fracture.dat, outcome.dat, fracture.mod,fracture.run
    Run the fracture.run file from AMPL. 
        It will generate two files: glow_max.txt (including maximum value of z) and glow_min.txt (including minimum value of z);
            
    Using an excel spreadsheet, calculate P_min from Z_max and P_max from Z_min (see the attached excel sheet)

### To reproduce results in Case Study 2:

    To solve maximization problem:

    Load all the files (DAT, MOD, RUN) from MAX folder in AMPL.

    Run the matched.run file from AMPL. 
        It will generate a file : Out_file_3stage (UB and LB for each N -total number of matched pair, and maximum value of LB at the end of each stage) 
        Record the maximum value of LB in the final stage as a final value for Z_max for each N
            
    To solve minimization problem:

    Load all the files (DAT, MOD, RUN) from MIN folder in AMPL.

    Run the matched.run file from AMPL. 
        It will generate a file : Out_min_file_3stage (UB and LB for each N -total number of matched pair, and minimum value of UB at the end of each stage)
        Record the minimum value of UB in the final stage as a final value for Z_min for each N

    Using an excel spreadsheet, calculate P_min from Z_max and P_max from Z_min (see the attached excel sheet)
