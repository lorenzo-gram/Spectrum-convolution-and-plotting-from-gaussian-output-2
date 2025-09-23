#!/bin/bash
module load python/3.8.12

python3 convolution_of_spectra.py --emin 1.5 --emax 3.1 --ewid 0.15 --grdid 300 --directories './*/*.log' --title './*/*.log' --transition 'abs' --experimental 'exp.dat'  

