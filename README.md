# OB_model
This is a biophysical model consisting of olfactory bulb mitral and granule cells. The granule cells are modeled as a dendritic and somatic compartment. The dendritic compartment forms dendrodendritic synapses with a subpopulation of the mitral cells. The somatic compartment receives input from a proportion of mitral cells and top down input modeled simplisticly. The mitral cell receives external input modeled as constant currents.

## Dependency
Code is designed to run in MATLAB R2018a. Some functions require `Curve Fitting Toolbox`. Can only run after MATLAB R2014b because of the `histcounts` function.

## How to use
- Download and unzip all the files
- The parameters are stored in `OB_params_GCE.txt` and can be updated.    

- To run it once and plot the power spectrum of simulated LFP oscillations, run the script `just_run.m`. The `input_file` is by default `OB_params_GCE.txt`. Make sure the directory is in MATLAB's path.     

- To run it with single/multiple-parameter sweep, run the script `run_with_sweep.m` and follow the detailed instruction inside the script.    

## Main functions
`InitNetwork_GCE.m`         -       initialize simulation parameters.     

`OB_network_GCE.m`          -       generate and initialize structs for neurons and network parameters.   

`NeuroActivity.m`           -       run the simulation, update the values for the neurons structs and generate a struct for synaptic current.   

`IandVLFP_GCE.m`            -       wrapper function for generating LFP from the outputs of `NeuroActivity.m`


## Author & Acknowledgement
This model is adapted from the one used by Osinski & Kay (2015). https://github.com/boleszek/GCExcitability   

The main function for network activity is broken into two files (`OB_network_GCE.m`, `NeuroActivity.m`) to fascilitate parameter sweep. They are modified to include a second compartment of the granule cell and top down input from the piriform cortex and cholinergic modulation. The numerical integration scheme is changed for better stability.

Besides these, with the exception of `InitNetwork_GCE.m`, `RasterPlot.m`, `subplot_tight.m`, `subplot1.m`, and `tightfig.m`, all the scripts and functions are written by Sam Zheng.
