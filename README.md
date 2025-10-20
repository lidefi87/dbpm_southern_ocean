# Dynamic Benthic Pelagic Model (DBPM) calibration - ISIMIP3A protocol
This repository contains all code necessary to process inputs used by DBPM. This repository has been redesigned to use both Python and R as part of the model workflow. Following protocol ISIMIP3A, this simulation uses inputs from GFDL-MOM6-COBALT2 at two horizontal resolutions: $0.25^{\circ}$ (original) and $1^{\circ}$ (coarsen).  
  
## Step 1. Processing DBPM climate inputs at a global scale
- Script [`01_processing_dbpm_global_inputs.ipynb`](scripts/01_processing_dbpm_global_inputs.ipynb) processes environmental data needed to force the DBPM model at a global scale. GFDL-MOM6-COBALT2 output files are transformed from `netCDF` to analysis ready `zarr` files. Files for `spinup` period are also created here.  
  
## Step 2. Processing DBPM climate inputs at a regional scale
- Script [`02_processing_dbpm_regional_inputs.ipynb`](scripts/02_processing_dbpm_regional_inputs.ipynb) uses `zarr` files produced in the previous step to extract data for an area of interest. In this notebook, we concentrate on the Southern Ocean, which was subdivided using three FAO Major Fishing Areas:
    - [FAO Major Fishing Area 48: Atlantic, Antarctic](https://www.fao.org/fishery/en/area/fao:48/en) referred to here as Weddell,  
    - [FAO Major Fishing Area 58: Indian Ocean, Antarctic And Southern](https://www.fao.org/fishery/en/area/fao:58/en) referred to here as East Antarctica,  
    - [FAO Major Fishing Area 88: Pacific, Antarctic](https://www.fao.org/fishery/en/area/fao:88/en) referred to here as West Antarctica.  

## Step 3. Processing DBPM fishing inputs at a regional scale
- Script [`03_processing_effort_fishing_inputs.R`](scripts/03_processing_effort_fishing_inputs.R) processes fishing catch and effort data for the area of interest. It also creates a single file including fishing and climate data, which has all variables needed to run DBPM within the boundaries of the area of interest.

## Step 4. Calculating fishing mortality parameters
- Script [`04_calculating_dbpm_fishing_params.R`](scripts/04_calculating_dbpm_fishing_params.R) does the following:  
    - Estimates fishing mortality parameters (catchability and selectivities for each functional group)  
    - Checks and adjusts the `search volume` parameter  
    - Creates and saves calibration plots in PDF format  
Plots created in this script can be used to visually inspect the fit of predicted catches against observed (reconstructed) catch data.
- Script [`04b_dbpm_runs_ccamlr_effort`](scripts/04b_dbpm_runs_ccamlr_effort.R) calculates fishing mortality parameters for DBPM runs using regional effort obtained from CCAMLR.

## Step 5. Setting up gridded inputs for spatial DBPM
- Script [`05_setup_gridded_DBPM.ipynb`](scripts/05_setup_gridded_DBPM.ipynb) processes all inputs necessary to run the spatial DBPM for the area and time period of interest.
- Script [`05a_setup_gridded_params_weddell.ipynb`](scripts/05a_setup_gridded_params_weddell.ipynb) applies a correction to inputs used to run DBPM in the Weddell Sea at $0.25^{\circ}$. Correction of inputs is needed for this region and resolution to deal with numerical instabilities.
- Script [`05b_setup_new_effort.ipnyb`](scripts/05b_setup_new_effort.ipynb) prepares data needed to run the gridded DBPM using CCAMLR effort data.  
  
## Step 6. Running DBPM spatial model  
- Script [`06_running_gridded_DBPM.ipynb`](scripts/06_running_gridded_DBPM.ipynb) uses inputs prepared in [step 5](scripts/05_setup_gridded_DBPM.ipynb) and runs the spatial DBPM. DBPM model outputs are stored for each timestep included in the input data.  

## Step 7. Calculating catches from gridded DBPM outputs 
- Script[`07_calculating_catches_gridded_DBPM.py`](scripts/07_calculating_catches_gridded_DBPM.py) calculates catches for benthic detritivores and pelagic predators from gridded DBPM outputs calculated in [step 6](scripts/06_running_gridded_DBPM.ipynb). Catch data is summarised per decade and maps created for the last decade of the spinup and the modelled period (1950 and 2010). Mean yearly catches are calculated for the area of interest from monthly catch estimates to create a time series.

## Step 8. Calculating biomass for predators and detritivores
- Script [`08_calculating_biomass.ipynb`](scripts/08_calculating_biomass.ipynb) produces calculates biomass for predators and detritivores from DBPM outputs produced in [step 6](scripts/06_running_gridded_DBPM.ipynb).  

## Step 9. Plotting outputs
- Script [`09_plotting_gridded_DBPM_outputs.ipynb`](scripts/09_plotting_gridded_DBPM_outputs.ipynb) creates plots and maps of biomass and catches calculated from DBPM outputs produced in [step 6](scripts/06_running_gridded_DBPM.ipynb).  
  
## Step 10. DBPM output evaluation
- Script [`10_evaluating_DBPM_outputs.R`](scripts/10_evaluating_DBPM_outputs.R) applies the observation range adjusted method to evaluate model performance as described in [Evans & Imran 2024](https://doi.org/10.1088/2515-7620/ad5ad8) prior to calculating a number of skill assessment metrics described in [Rynne et al 2025](https://doi.org/10.1029/2024EF004868).  

## Step 11. Comparing effort datasets
- Script [`11_effort_plots`](scripts/11_effort_plots.R) compares effort data as prescribed by the FishMIP Protocol 3a and CCAMLR data.  

## Step 12. Plots for publication
- Script [`12_supporting_plots`](scripts/12_supporting_plots) produces plots to support contents of manuscript describing exploration on DBPM outputs within the Southern Ocean.  

The `useful_functions` scripts in Python and R contain functions that have been developed to run DBPM.
  
# Running this repository
The scripts in this repository were developed in NCI's Gadi, so the easiest way to run these script is to clone this repository to Gadi. However, before you can do this, you will need an NCI account, which are only available for researchers with an email address from an Australian institution. Further down in this document, we include information about how to create an NCI account if you do not have one already. Remember, you must have an email address for an Australian institution to create an NCI account.  
  
You can also run these scripts in your own computer or a different server, but you will need need access to the forcing data (i.e., GFDL-MOM6-COBALT2 outputs and fishing data) to run them. We include information about how to access these data.  
  
## Getting an NCI account
1. [Create an NCI user account](https://access-hive.org.au/getting_started/first_steps/#create-an-nci-user-account)  
      * You should use your Australian institution’s email account when signing up  
      * When it asks for a project to join:  
        * If possible, contact the NCI scheme manager at your institution to find out what NCI project you should use to sign up for your NCI account. This account will provide you with computing power.    
2. [Join relevant NCI projects](https://access-hive.org.au/getting_started/first_steps/#join-relevant-nci-projects)
      * Request to join the following NCI projects:  
        * vf71 - for access to GFDL-MOM6-COBALT2 outputs in analysis ready data format 
        * xp65 - for the Python conda environment   
      * Note that it can take a few business day get approved as a project member
3. [Verify that you can log into NCI’s Gadi](https://access-hive.org.au/getting_started/first_steps/#login-to-gadi)  
      * Note that it usually takes more than 30 minutes for your account to be created  
      * You are also welcome to follow along with the [instructions to set up ssh keys](https://access-hive.org.au/getting_started/first_steps/#automatic-login), but this is optional.  

## Accessing forcing data 
### Ocean outputs from GFDL-MOM6-COBALT2
The environmental data comes from GFDL-MOM6-COBALT2, which is available at two horizontal resolutions: $0.25^{\circ}$ (original model outputs) and ($1^{\circ}$, coarsen from original outputs). The original GFDL-MOM6-COBALT2 outputs can be downloaded from the [Inter-Sectoral Impact Model Intercomparison Project (ISIMIP) Data Portal](https://data.isimip.org/search/tree/ISIMIP3a/InputData/climate/ocean/gfdl-mom6-cobalt2/) as `netCDF` files. However, you can also access GFDL-MOM6-COBALT2 outputs as `zarr` files from project `vf71` at the [National Computational Infrastructure (NCI)](https://nci.org.au/).  
  
### Fishing catch data
The fishing catch data came from three sources:  
1. 'ISIMIP3a reconstructed fishing activity data (v1.0)' ([Novaglio et al. 2024](https://data.isimip.org/10.48364/ISIMIP.240282))  
2. 'CCAMLR Statistical Bulletin, Vol. 36' ([CCAMLR 2024](https://www.ccamlr.org/en/publications/statistical-bulletin))  
3. 'Sea Around Us catch reconstructions' ([Pauly et al. 2020](https://www.seaaroundus.org/data/#/fao))

### CCAMLR effort data
The regional effort data used in comparisons came from the CCAMLR Statistical Bulletin, Vol. 36. This data is publicly available and can be downloaded from their [website](https://www.ccamlr.org/en/publications/statistical-bulletin).  
  
Copies of these datasets are also available under project `vf71` at the [National Computational Infrastructure (NCI)](https://nci.org.au/).  

