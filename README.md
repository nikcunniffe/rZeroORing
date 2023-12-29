# rZeroORing
Code for paper by van den Bosch, Helps and Cunniffe on using the O-ring statistic to calculate R0 for spatially structured epidemic models

0/ Compile EpidemicSim.exe from EpidemicSim.c and mt19937ar.c
1/ Create directory to do the runs
2/ Copy the following files to directory created in step 0
	- EpidemicSim.cfg
	- EpidemicSim.exe
	- create_LS.R
	- rZero_From_Sims.R
	- rZero_Function.R
3/ Create the following subdirectories of directory created in step 0
	- Inputs
	- Outputs
4/ Start R and change to the directory created in step 0 using setwd()
5/ Run create_LS.R 
	- Options for landscape generation are in the R file
	- Running it will fill up Inputs subdirectory
6/ Run EpidemicSim.exe on command line
	- Options for epidemics are in the EpidemicSim.cfg files
	- will fill up Outputs subdirectory
7/ Run rZero_From_Sims.R
	- will print estimated and calculated rZero to the screen
