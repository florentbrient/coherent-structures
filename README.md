# Coherent Structures in atmopsheric boundary layers in Large Eddy Simulations
__This project aims to characterize coherent structures in LES.__  
Git project creation in 2019  

The goal of ...

The routines have been used in a peer-reviewed article :

Brient, Couvreux, Rio, Honnert 2023, *Coherent subsiding structures in large eddy simulations of atmospheric boundary layers*, Quarterly Journal of the Royal Meteorological Society (DOI: TO ADD). See preprint article on [EarthAriXiv](https://eartharxiv.org/repository/view/5713/) 

Authors : F. Brient, F. Couvreux, C. Rio, R. Honnert

## Description

### Main directory
The main repository contains this __README.md__ file and XX folders: __rename_netcdf__ and __scripts__.

The main script is __exec_step.exe__, run as *./exec_step.exe 1 1 1 1 1*

It calls the following routines:

### 1. Rename variables
The first call runs __script_rename_dim__ on the *rename_netcdf/* folder

It modifies the outputs of the MESO-NH to remove unnecesseray variables and dimensions.

### 2. Object selection
The object selection use a collaborative project aimed develop generic tools to identify and characterize objects in high resolution model output fields. It is available here: [https://gitlab.com/tropics/objects](https://gitlab.com/tropics/objects)

### 3. Characterize objects
The routine aims to characterize simulations by providing figures of *x-z* and *x-y* cross sections, and domain-averaged vertical profiles.

Some options can be chosen on the __plot_LES_3D.py__:
- Activate/deactivate objects
- Mean flux or characteristics
- Threshold of objects (*thrs*)
- Define inversion height (*idxzi*)
- Subdomain for cross sections
- Filtering cloudiness
- Averaging over a layer

### 4. Calculate statistics of objects









