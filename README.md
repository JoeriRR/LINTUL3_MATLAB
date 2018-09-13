# LINTUL3_MATLAB
LINTUL3 crop model in MATLAB obtained from Wagingen University. The FORTRAN code inside the FST-simulator has been transferred to MATLAB code.
## Installation
Just clone this repository, and change the path 'weather_file' path variable to the correct location to run LINTUL3.m. to Run LINTUL3_res_clean.m make sure that res.dat is in the same folder as the script.
## Overview of files
This chapter describes the files.
### LINTUL3.m
Script to run the LINTUL3 crop model inside the MATLAB environment.
### LINTUL3_res.m
Script to compare the results from MATLAB and the FST-simulator. 
### monthday.m
Function to obtain month and day.
### NLD6.987
Weatherfile 1987, Amersfoort, Netherlands.
### importfile.m
Function to read res.dat.
### res.dat
Data obtained from the FST-simulator.
### Read_Weatherfile.m
Function to read the weatherfile NLD6.987.
### Z.mat 
File which has all the same variables as res.dat but from MATLAB instead of FST.
