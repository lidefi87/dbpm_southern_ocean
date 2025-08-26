## Running DBPM with CCAMLR effort 
# Is effort data responsible for the differences between observed and estimated
# catches?


# Loading libraries -------------------------------------------------------
source("scripts/useful_functions.R")
library(arrow)
library(dplyr)
library(ggplot2)
library(cowplot)
library(stringr)
library(grid)
library(gridExtra)
library(jsonlite)
library(zoo)
library(purrr)


# Preparing CCAMLR effort to force model ----------------------------------

# Load effort data from CCAMLR
eff_ccamlr_int <- read_parquet("data/effort_per_year_by_area.parquet") |> 
  # Keep effort data up to 2010
  filter(year <= 2010) |> 
  # Ensure there is one row for every year between first and last effort record
  group_by(asd_area_code) |> 
  complete(year = seq(min(year), max(year))) |>
  # For any missing years, interpolate effort data 
  mutate(effort = na.approx(Effort_kW_Days_yr),
         region = str_c("FAO ", asd_area_code)) |> 
  ungroup() |> 
  select(region, year, effort)

# Load original inputs - including original effort forcing from Novaglio
eff <- list.files("/g/data/vf71/la6889/dbpm_inputs/", pattern = "^[e|w]",
                  full.names = T) |> 
  map_chr(\(x) list.files(file.path(x, "monthly_weighted/025deg"), 
                          pattern = "dbpm_clim-fish-inputs_fao", 
                          full.names = T)) |> 
  map(\(x) read_parquet(x)) |> 
  bind_rows()

# Replace original effort data with CCAMLR effort data when available
new_eff <- eff |> 
  left_join(eff_ccamlr_int, by = c("region", "year")) |>
  # Calculate relative effort
  group_by(region) |> 
  mutate(eff_area = effort/tot_area_m2,
         eff_rel_area = eff_area/max(eff_area, na.rm = T),
         nom_active_area_m2_relative = case_when(is.na(eff_rel_area) ~ 
                                                   nom_active_area_m2_relative,
                                                 T ~ eff_rel_area)) 

# Plotting data to ensure effort was correctly updated
new_eff |> 
  ggplot()+
  geom_line(aes(year, eff_rel_area, color = as.factor(region)))+
  geom_line(aes(year, nom_active_area_m2_relative, color = as.factor(region)), 
            linetype = "dashed")+
  facet_grid(region~.)

# Save new effort data 
new_eff |> 
  write_parquet("data/new_effort_ccamlr_added.parquet")


# Run non-spatial DBPM ----------------------------------------------------
# Define base directory
base_dir <- "/g/data/vf71/la6889/dbpm_inputs"

# Define variables for region of interest
region_int <- "fao-48"
region_name <- "weddell"
reg_name_plot <- "Atlantic Sector"
res <-  "1deg"
base_folder <- file.path(base_dir, region_name)
dbpm_out_folder <- file.path(base_folder, "run_nonspatial")

# Subset input data with CCAMLR effort for region of interest
dbpm_inputs <- new_eff |> 
  filter(region == str_replace(str_to_upper(region_int), "-", " "))

# Searching best fishing parameters values for area of interest -----------
# Path to folder where results will be stored
results_folder <- file.path(base_folder, "fishing_params", 
                            paste0("best_fish_vals_", region_int, "_", res))

# Find fishing parameters used for 
fishing_params <- list.files(results_folder, "best-fishing-parameters",
                             full.names = T) |> 
  read_parquet() |> 
  arrange(rmse) |> 
  slice(1)
  
# Create non-spatial parameters
params <- sizeparam(dbpm_inputs, fishing_params, xmin_consumer_u = -3,
                    xmin_consumer_v = -3, tstepspryr = 12)

# Save parameters - ensure 10 decimal places are saved
params |> 
  write_json(file.path(results_folder, 
                       paste0("dbpm_size_params_ccamlr_effort_", region_int, 
                              ".json")), digits = 10)

# Run non-spatial DBPM
non_spatial_run <- run_model(fishing_params, dbpm_inputs)

# Save results
non_spatial_run |>
  write_parquet(file.path(dbpm_out_folder,
                          paste0("dbpm_nonspatial_ccamlr_effort_",region_int,
                                 "_1841-2010.parquet")))


# Comparing catches from new and old effort -------------------------------
# Original catch estimates
catch_dbpm_nonspatial <- list.files(file.path(base_folder, 
                                              "run_nonspatial/1deg"), 
                                    full.names = T) |> 
  read_parquet() |> 
  filter(year >= 1961) |> 
  select(c(time, year, total_catch)) |> 
  group_by(year) |> 
  summarise(`effort_Novaglio et al` = mean(total_catch, na.rm = T))

# Observed catches
catch_obs <- list.files(file.path(base_folder, "monthly_weighted/1deg"),
                        pattern = "^dbpm_clim", full.names = T) |> 
  read_parquet() |> 
  filter(year >= 1961) |> 
  group_by(year) |> 
  summarise(obs_ccamlr = mean(catch_ccamlr*1e6, na.rm = T),
            obs_pauly = mean(catch_pauly*1e6, na.rm = T),
            obs_watson = mean(catch_tonnes_area_m2*1e6, na.rm = T)) |> 
  rowwise() |>
  mutate(min_catch_density = min(obs_ccamlr, obs_pauly, obs_watson, na.rm = T),
         max_catch_density = max(obs_ccamlr, obs_pauly, obs_watson, 
                                 na.rm = T)) |> 
  select(year, min_catch_density, max_catch_density)

# New catch estimates - CCAMLR effort
catch_dbpm_ccamlr <- non_spatial_run |> 
  filter(year >= 1961) |> 
  group_by(region, year) |> 
  summarise(effort_CCAMLR = mean(total_catch, na.rm = T))

# Merge all catches into single data frame
catch_all <- catch_dbpm_ccamlr |> 
  full_join(catch_dbpm_nonspatial, by = "year") |> 
  full_join(catch_obs, by = "year") |> 
  pivot_longer(cols = starts_with("effort_"), names_to = "source", 
               names_prefix = "effort_", values_to = "vals")


# Plots comparing catch estimates -----------------------------------------
# Define variables for region of interest
region_int <- "fao-88"
region_name <- "west_antarctica"
reg_name_plot <- "Pacific Sector"
  
catch_all_reg <- read_parquet(
  "data/catch_comp_ccamlr_novaglio_all_regions.parquet") |> 
  mutate(min_catch_density = ifelse(is.infinite(min_catch_density), NA,
                                    min_catch_density),
         max_catch_density = ifelse(is.infinite(max_catch_density), NA,
                                    max_catch_density)) |> 
  group_by(region) |> 
  group_split()

catch_88 <- catch_all_reg[[3]] |> 
  ggplot(aes(x = year))+
  geom_ribbon(aes(ymin = min_catch_density, ymax = max_catch_density), 
              fill = "#e3dded")+
  geom_line(aes(y = vals, linetype = source), color = "#fdcc8a", 
            linewidth = 0.9)+
  labs(linetype = "Fishing effort source", 
       y = ~paste("Mean annual catch per unit area (g*", m^-2, ")"),
       subtitle = paste0(str_to_title(reg_name_plot), " (", 
                         str_replace_all(str_to_upper(region_int), "-", " "),
                         ")"))+
  theme_bw()+
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(),
        axis.title.y = element_text(family = "sans", size = 14),
        axis.text = element_text(family = "sans", size = 12),
        legend.title = element_text(face = "bold", size = 12),
        legend.text = element_text(family = "sans", size = 12),
        plot.subtitle = element_text(family = "sans", size = 12, hjust = 0.5))+
  theme(legend.position = "none", axis.title.y = element_blank())

leg <- get_legend(catch_48)

catch_48 <- catch_48+
  theme(legend.position = "none", axis.title.y = element_blank())

ytitle <- textGrob(~paste("Mean annual catch per unit area (g*", m^-2, ")"),
                   rot = 90)

fig <- grid.arrange(arrangeGrob(plot_grid(plot_grid(catch_48, catch_58, 
                                                    ncol = 2), 
                                          plot_grid(catch_88, leg), nrow = 2),
                                left = ytitle))

# legend.byrow = T, legend.direction = "horizontal", 
# legend.position = "bottom",
# 
# plot.margin = margin(5, 15, 5, 5, "pt"),





