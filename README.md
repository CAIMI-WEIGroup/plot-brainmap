# plot-brainmap
Plot brain maps based on freesurfer fsaverage template

# Examples
Create a random vector to plot, base on the DK114 atlas ('lausanne120')

val = randn(114, 1)
plotBrainRegions(val, 'lausanne120');

# References
cbrewer: 
Scott Lowe (2022). cbrewer2 (https://github.com/scottclowe/cbrewer2), GitHub. Retrieved October 28, 2022.
real2rgb: 
Oliver Woodford (2022). real2rgb & colormaps (https://www.mathworks.com/matlabcentral/fileexchange/23342-real2rgb-colormaps), MATLAB Central File Exchange. Retrieved October 28, 2022.