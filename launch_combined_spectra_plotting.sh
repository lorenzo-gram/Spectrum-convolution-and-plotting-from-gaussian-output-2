#!/bin/bash
module load python/3.8.12

python3 combined_spectra_plotting.py --emin 1.5 --emax 3.1 --ewid 0.15 --grdid 300 --directories_gaussian './*/*.log' --directories_orca './*/*.out' --title_gaussian './*/*.log' --title_orca './*/*.out' --transition 'abs' --experimental 'exp.dat'  

