#!/usr/bin/env python3

# Calculating catches from gridded DBPM outputs
# Authors: Denisse Fierro Arcos
# Date last update: 2025-08-06

#Loading libraries
import os
from glob import glob
import xarray as xr
import useful_functions as uf
import json

#Name of region and model resolution
reg_number = 'fao-58'
reg_name = 'east_antarctica'
model_res = '025deg'
# Choose between '_simask' and '_ccamlr_eff'
runs = '_simask'
# Optional - Subset biomass for predators and detritivores
# If different from the minimum body size given in the gridded_params variable
# If not needed, set to "None", otherwise select size class
sub_sc = 10

#Defining input and output folders
base_folder = f'/g/data/vf71/la6889/dbpm_inputs/{reg_name}'
gridded_inputs = os.path.join(base_folder, 'gridded_params', model_res)
if runs == '_simask':
    gridded_outputs = os.path.join(base_folder, 'run_fishing_seaicemask', 
                                   model_res)
elif runs == '_ccamlr_eff':
    gridded_outputs = os.path.join(base_folder, 'run_fishing_new_eff', 
                                   model_res)
else:
    gridded_outputs = os.path.join(base_folder, 'run_fishing', 
                                   model_res)
outputs_folder = os.path.join(base_folder, 'gridded_dbpm_outputs', model_res)

#Ensure outputs folder exists
os.makedirs(outputs_folder, exist_ok = True)

## Loading gridded parameters and gridded inputs ---
#Gridded DBPM parameters
gridded_params = json.load(open(
    os.path.join(gridded_inputs, f'dbpm_gridded_size_params_{reg_number}_python.json')))

#Size class bins
log10_size_bins_mat = xr.open_zarr('../outputs/log10_size_bins_matrix.zarr/')['size_bins']
size_bin_vals = (10**log10_size_bins_mat)

#Area - to be used in weighted mean calculation
area = xr.open_zarr(glob(os.path.join(base_folder, 'gridded', model_res, 
                                      '*areacello*'))[0])['cellareao']
area = area.fillna(0)

## Calculate biomass for predators and detritivores ---
#Predators
var_pred = 'predators'
pred_raw = (xr.open_mfdataset(glob(
    os.path.join(gridded_outputs, f'{var_pred}*[0-9][0-9][0-9].nc')))[var_pred])

#Detritivores
var_det = 'detritivores'
det_raw = (xr.open_mfdataset(glob(
    os.path.join(gridded_outputs, f'{var_det}*[0-9][0-9][0-9].nc')))[var_det])

pred = ((pred_raw*size_bin_vals*
        gridded_params['log_size_increase']).
    isel(size_class = slice(gridded_params['ind_min_pred_size'],
                            None)))
pred.name = 'pred_biomass'
pred['size_class'] = 10**pred.size_class.values
pred = pred.assign_attrs({'short_name': 'pred_biomass',
                          'long_name': 'Predator biomass per size class (in g)',
                          'units': 'g m-2'})

det = ((det_raw*size_bin_vals*
        gridded_params['log_size_increase']).
    isel(size_class = slice(gridded_params['ind_min_detritivore_size'],
                            None)))
det.name = 'detriti_biomass'
det['size_class'] = 10**det.size_class.values
det = det.assign_attrs({'short_name': 'detriti_biomass',
                        'long_name': 'Detritivore biomass per size class (in g)',
                        'units': 'g m-2'})
## Save outputs ---
pred.to_zarr(
    os.path.join(outputs_folder,
                 f'pred_biomass_size_group{runs}_{reg_number}_1841-2010.zarr'),
    consolidated = True, mode = 'w')

det.to_zarr(
    os.path.join(outputs_folder,
                 f'detriti_biomass_size_group{runs}_{reg_number}_1841-2010.zarr'),
    consolidated = True, mode = 'w')

## Optional: Subset biomass for predators and detritivores ---
# If different from the minimum body size given in the gridded_params variable
if sub_sc is not None:
    pred = pred.sel(size_class = slice(sub_sc, None))
    det = det.sel(size_class = slice(sub_sc, None))
    
## Calculate mean group biomass per area per time step ---
# Predators
pred_tot = pred.sum('size_class')
pred_tot_w = pred_tot.weighted(area)
pred_tot_w_mean = pred_tot_w.mean(('lat', 'lon')).load()
pred_tot_w_mean.name = 'tot_pred_bio'

# Detritivores
det_tot = det.sum('size_class')
det_tot_w = det_tot.weighted(area)
det_tot_w_mean = det_tot_w.mean(('lat', 'lon')).load()
det_tot_w_mean.name = 'tot_det_bio'

## Calculate mean group exploitable biomass per area per time step ---
#Predators
pred_exp = ((pred_raw*size_bin_vals*
             gridded_params['log_size_increase']).
    isel(size_class = slice(gridded_params['ind_min_fish_pred'],
                            None))).sum('size_class')
pred_exp_w = pred_exp.weighted(area)
pred_exp_w_mean = pred_exp_w.mean(('lat', 'lon')).load()
pred_exp_w_mean.name = 'tot_expl_pred_bio'

#Detritivores
det_exp = ((det_raw*size_bin_vals*
             gridded_params['log_size_increase']).
    isel(size_class = slice(gridded_params['ind_min_fish_det'],
                            None))).sum('size_class')
det_exp_w = det_exp.weighted(area)
det_exp_w_mean = det_exp_w.mean(('lat', 'lon')).load()
det_exp_w_mean.name = 'tot_expl_det_bio'

## Saving results as tables ---
bio_pred_det = (xr.merge([pred_exp_w_mean, pred_tot_w_mean, 
                          det_exp_w_mean, det_tot_w_mean]).to_dataframe().
    reset_index())
bio_pred_det['region'] = reg_number.upper().replace('-', ' ')
bio_pred_det['resolution'] = model_res

yr_start = bio_pred_det.time.dt.year.min()
bio_pred_det.to_parquet(
    os.path.join(outputs_folder,
                 # f'mean_pred_det_bio_{reg_number}_1841-2010.parquet'),
                 f'mean_pred_det_10g-1t-bio_exploit-bio{runs}_{reg_number}_{yr_start}-2010.parquet'),
    index = False)
