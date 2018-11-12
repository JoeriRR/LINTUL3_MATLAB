# LINTUL3_MATLAB
LINTUL3 crop model in MATLAB obtained from Wagingen University. The FORTRAN code inside the FST-simulator has been transferred to MATLAB code (LINTUL3.m/LINTUL3_onestage.m). This MATLAB code is compared with the original FORTRAN code in LINTUL3_res.m, which uses the file 'res.dat' as input. The weatherfiles for Amersfoort that are used for input can be found in the ICWEATHER directory. Results from the DHMARA optimization can be found in the Results directory. Figures which are saved in the GUI are saved in the directory Figures.
## Installation
Just clone this repository and it should work. 
## Files overview
Each file is described below.
### LINTUL3 scripts
The scripts LINTUL3.m and LINTUL_onestage simulate the LINTUL3 model for the whole season and daily. The script LINTUL_res.m compares the MATLAB LINTUL3 model with the FST-simulator for the crop spring-wheat.
### Optimization simulation scripts
The scripts that start with DHMARA simulate the optimization problem with inventory constraints (refill) with spatial constraints (spatial) with virtual field (virt_field) and the original optimization (with_fert). 
### Read scripts
The scripts importfile.m and Read_Weatherfile read the 'res.dat' file and weatherdata and loads them as a matrix.
### Plot scripts
The scripts monthday.m and plotyear.m are used for plots.
### Test scripts
The scripts test_inventory_control.m and test_virtual_field are scripts to validate the inventory and virtual field constraints.
