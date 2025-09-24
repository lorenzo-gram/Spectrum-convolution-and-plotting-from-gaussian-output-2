# Spectrum-convolution-and-plotting-from-gaussian-output

This python code is useful to plot the result of one or various tddft calculations performed by Gaussian and, if needed, compare it with an experimental spectrum.  
It works with either one output file or more output files. In both cases, for each output file indicated, the script will plot a gaussian function around each excitation (whose intensity is proportional to the oscillator strength in case of absorption and proportional to the [oscillator strength]*[emission wavelength]^2 in the case of the emission) energy and sum up all the gaussian functions.  
If more than one output files are provided, it will plot the different spectra (one for each output file) together (each one with a different color).
It is also possible to plot the experimental spectrum for comparison: the wavelength and intensity must be in a .dat file in two separate columns (no title on top of the columns).

## How to launch the code?

The code can be launched directly on the command line or, for the laziest like me, it is possible to write a simple .sh script so that every time you need to use the script you can simply modify the .sh :) An example is given here "launch_convolution_of_spectra_separate.sh".  

The code is launched with:  

python3 convolution_of_spectra.py --emin 1.5 --emax 4.5 --ewid 0.15 --grdid 300 --directories './md*/OUTPUT/QM_data/qmALL.log' --title './md*' --transition 'abs' --nstates 5

where:
* --emin specifies the minimum energy to be plotted (in eV)
* --emax specifies the maximum energy to be plotted (in eV)
* --ewid 0.15 specifies the standard deviation of each gaussian function (in eV)
* --grdid specifies the number of points in the plot of the spectrum
* --directories specifies (in python string format) the path to the .log file(s)
* --title specifies the path where you find the directories that will give the name to each different gaussian calculation (in this case md001, md002, md003,...)
* --transition specifies if we are simulating absorption ('abs') or emission ('emi') spectrum
* --nstates specifies the number of transitions to plot and put in the legend
