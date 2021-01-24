
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