# Replication Code for the Paper: A robust approach to quantifying uncertainty in matching problems of causal inference

Link to the paper: https://arxiv.org/abs/1812.02227

*R packages required to run the code:* 
- reshape2, ggplot2, ggjoy, latex2exp, plyr

*Python packages required to run the code:*
- numpy, scipy, pandas (for Figures 2 and 10)
- numpy, scipy, matplotlib, tabulate (for Figure 9)

*Software required: AMPL, Excel*

### To reproduce results in Figures 2 and 10 (appendix) 
    run: simulations/simulate_and_plot.R

### To reproduce results in Figure 9 (appendix) 
    run: distribution_comparison.py

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
