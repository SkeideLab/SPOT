# 01_dataprep

Data preparation consists of alignment of native surfaces to group template, warping of retinotopy atlases via the group template into individual surface space, projecting functional data to the native surfaces, and the simulation of alternative surface data based on a model of local spatial correlations. 

## Why alignment and warping?

The main analysis is being conducted in native surface space of each participant to avoid warping and interpolating the data as much as possible. However, the retinotopy atlases are in fsaverage space. To get the template to each individual surface space, the surfaces are registered to each other via the dHCP group surface. First, the dHCP group surface dhcpSym is registered to fsaverage via the existing registrations of dhcpSym to HCP's fs_lr and of HCP's fs_lr to fsaverage. Then, for each participant, their native surface is registered to the dHCP group space. The native -> dhcpSym and dhcpSym -> fsaverage registrations are then combined to make one native -> fsaverage registration. 

## Directory structure

`alignment`: contains code and settings files for alignment of native surfaces to a group surface. Adapted from [dhcp template alignment](https://github.com/ecr05/dHCP_template_alignment).

`retinotopy`: contains code to produce retinotopy gifti files from mgz, to warp the templates to native space, to partition the visual areas, the retinotopy templates themselves, and registrations between standard spaces necessary to do the warping.

`sanitycheck`: code for a quick default mode network analysis on the functional surface data to check their integrity.

`simulation_alternativemodel`: code to produce simulated surface data based on a model of random data and local spatial correlations.

`surfaceprojection`: code to project functional data to the surface.

`run_surfaceprep.sh`: code to run everything in this directory and prepare all the data. 