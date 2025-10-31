#!/usr/bin/env python3

# Processing DBPM outputs before creating maps
# Authors: Denisse Fierro Arcos
# Date last update: 2025-10-30

#Loading libraries
import os
from glob import glob
import xarray as xr
import re
import numpy as np
import json

## Model resolution and run identifier ----
model_res = '1deg'
# Can be either '_simask' or blank ('')
runs = '_simask'

## Get list of relevant files
base_folder = '/g/data/vf71/la6889/dbpm_inputs'
area_files = glob(os.path.join(base_folder, f'*/gridded/{model_res}/*areacello*'))
si_files = glob(os.path.join(base_folder, f'*/gridded/{model_res}/*obsclim_siconc*month*'))
catch_files = glob(
    os.path.join(base_folder, f'*/gridded_dbpm_outputs/{model_res}/total_catches{runs}_{model_res}_*'))
bathy_files = glob(os.path.join(base_folder, f'*/gridded/{model_res}/*obsclim_deptho_*'))
pred_files = glob(os.path.join(base_folder, 
                               f'*/gridded_dbpm_outputs/{model_res}/pred_*group{runs}_fao*1841*'))
det_files = glob(os.path.join(base_folder, 
                               f'*/gridded_dbpm_outputs/{model_res}/detriti_*group{runs}_fao*1841*'))
fishparam_files = glob(os.path.join(base_folder, 
                                    f'*/gridded_params/{model_res}/dbpm_gridded_size_params_*_python.json'))

#Size class bins
log10_size_bins_mat = xr.open_zarr('../outputs/log10_size_bins_matrix.zarr/')['size_bins']
size_bin_vals = 10**log10_size_bins_mat

# Create empty dictionary to store catch, depth and sea ice extent results for all sectors
catch_all = {}
catch_month_all = {}
sie_all = {}
sie_month_all = {}
bathy_all = {}
tot_exp_all = {}

# Looping through all sectors
for i, af in enumerate(area_files):
    # Identify sector from file name
    [region] = re.findall('fao-[0-9]{2}', af)
    # Create land mask from area files
    mask = xr.open_zarr(af)['cellareao']
    mask = xr.where(np.isfinite(mask), 1, np.nan)

    # Identify sea ice extent from sea ice concentration data
    si = xr.open_zarr(si_files[i])['siconc']
    #Calculate montly sea ice mean
    si_mean = si.mean('time')*mask
    # Convert to sea ice extent
    sie_mean = xr.where(si_mean >= 15, 1, 0)*mask
    # Adding sea ice extent to dictionary
    sie_all[region] = sie_mean

    #Calculating monthly sea ice extent
    si_month = si.groupby('time.month').mean('time')*mask
    sie_month = xr.where(si_month >= 15, 1, 0)*mask
    sie_month_all[region] = sie_month

    # Load estimated catches
    total_catch = xr.open_zarr(catch_files[i])['total_catches'].sel(time = slice('1961', None))
    # Calculate mean catches for entire model period (1961-2010)
    mean_catch = total_catch.mean('time')*mask
    mean_catch.name = 'mean_fishing_catches'
    # Adding catches to dictionary
    catch_all[region] = mean_catch

    #Calculating monthly mean catches
    mean_month_catch = total_catch.groupby('time.month').mean('time')*mask
    mean_month_catch.name = 'mean_monthly_fishing_catches'
    catch_month_all[region] = mean_month_catch

    # Load depth
    bathy_all[region] = xr.open_zarr(bathy_files[i])['deptho']

    # Load predators
    predators = xr.open_zarr(pred_files[i])['pred_biomass'].sel(time = slice('1961', None))
    
    #Loading detritivores
    detritivores = xr.open_zarr(det_files[i])['detriti_biomass'].sel(time = slice('1961', None))
    
    #Model parameters
    gridded_params = json.load(open(fishparam_files[i]))
    
    pred_exp = (predators.
        sel(size_class = slice(size_bin_vals[gridded_params['ind_min_fish_pred']], None)).
        sum('size_class'))
    
    det_exp = (detritivores.
        sel(size_class = slice(size_bin_vals[gridded_params['ind_min_fish_det']], None)).
        sum('size_class'))

    #Calculate total exploitable biomass
    tot_exp = (pred_exp+det_exp).mean('time')*mask
    tot_exp.name = 'total_exploitable_biomass'
    tot_exp_all[region] = tot_exp
    

# Transforming dictionaries to datasets
#Catches
catch_all = xr.Dataset(catch_all)
catch_all = catch_all.chunk({'lat': -1, 'lon': -1})
#Catches - monthly
catch_month_all = xr.Dataset(catch_month_all)
catch_month_all = catch_month_all.chunk({'lat': -1, 'lon': -1})

#Sea ice extent
sie_all = xr.Dataset(sie_all)
sie_all = sie_all.chunk({'lat': -1, 'lon': -1})
#Sea ice extent - monthly
sie_month_all = xr.Dataset(sie_month_all)
sie_month_all = sie_month_all.chunk({'lat': -1, 'lon': -1})

#Bathymetry
bathy_all = xr.Dataset(bathy_all)
bathy_all = bathy_all.chunk({'lat': -1, 'lon': -1})
#Exploitable biomass
tot_exp_all = xr.Dataset(tot_exp_all)
tot_exp_all = tot_exp_all.chunk({'lat': -1, 'lon': -1})

# Saving results
#Catches
catch_all.to_zarr(
    os.path.join('../data', f'mean_catch_all_regions{runs}_{model_res}_1961-2010.zarr'),
    consolidated = True, mode = 'w')
#Catches - monthly
catch_month_all.to_zarr(
    os.path.join('../data', f'monthly_mean_catch_all_regions{runs}_{model_res}_1961-2010.zarr'),
    consolidated = True, mode = 'w')

#Sea ice extent
sie_all.to_zarr(
    os.path.join('../data', f'mean_sea-ice-edge_all_regions_{model_res}_1961-2010.zarr'),
    consolidated = True, mode = 'w')
#Sea ice extent - monthly
sie_month_all.to_zarr(
    os.path.join('../data', f'monthly_mean_sea-ice-edge_all_regions_{model_res}_1961-2010.zarr'),
    consolidated = True, mode = 'w')

#Bathymetry
bathy_all.to_zarr(os.path.join('../data', f'depth_all_regions_{model_res}_fixed.zarr'),
                  consolidated = True, mode = 'w')
#Exploitable biomass
tot_exp_all.to_zarr(
    os.path.join('../data', f'mean_tot-exp-bio_all_regions{runs}_{model_res}_1961-2010.zarr'),
    consolidated = True, mode = 'w')
