
# Loading libraries -------------------------------------------------------
source("scripts/useful_functions.R")
library(ggplot2)
library(dplyr)
library(arrow)
library(jsonlite)
library(lubridate)
library(stringr)
library(purrr)
library(cowplot)


# Non-spatial exploited biomass -------------------------------------------
# Creating empty data frame to store exploited biomass results
exp_bio_data <- data.frame()
base_folder <- "/g/data/vf71/la6889/dbpm_inputs"

# Define variables identifying regions
region_int <- "fao-48"
region_name <- "weddell"

# Define location of non-spatial results
results_folder <- file.path(base_folder, region_name, "fishing_params", 
                            paste0("best_fish_vals_", region_int, "_1deg"))

# Load DBPM inputs
dbpm_inputs <- file.path(base_folder, region_name, "monthly_weighted/1deg",
                         paste0("dbpm_clim-fish-inputs_", region_int, 
                                "_1841-2010.parquet")) |> 
  read_parquet()

# Load all parameters used to run DBPM
params <- read_json(file.path(results_folder, 
                              paste0("dbpm_size_params_", region_int, ".json")), 
                    simplifyVector = T)

# Load fishing parameters used to run DBPM
fishing_params <- read_parquet(
  file.path(results_folder, 
            paste0("best-fishing-parameters_", region_int, 
                   "_searchvol_estimated_numb-iter_1000.parquet"))) |> 
  arrange(rmse) |> 
  slice(1)

# Run non-spatial DBPM
non_spatial_run <- run_model(fishing_params, dbpm_inputs, withinput = F)

# Calculate exploited predator and detritivore biomass
size_bin_vals <- (10^params$log10_size_bins)
pred_ind <- params$ind_min_fish_pred:length(size_bin_vals)
pred_exp <- colSums((non_spatial_run$predators*size_bin_vals*
                       params$log_size_increase)[pred_ind,])
det_ind <- params$ind_min_fish_det:length(size_bin_vals)
det_exp <- colSums((non_spatial_run$detritivores*size_bin_vals*
                      params$log_size_increase)[det_ind,])

# Create data frame with results
exp_bio <- data.frame(year = dbpm_inputs$year, region = fishing_params$region,
                      resolution = "non-spatial", pred_exp = pred_exp, 
                      det_exp = det_exp) |> 
  # Calculate mean annual biomass
  group_by(year, region, resolution) |> 
  summarise(across(ends_with("_exp"), list(mean = mean))) |> 
  # Add all exploited biomass and add region and resolution information
  mutate(mean_expl_bio = pred_exp_mean+det_exp_mean)

# Merging data per region into main data frame
exp_bio_data <- exp_bio_data |> 
  bind_rows(exp_bio)

rm(region_int, region_name, results_folder, dbpm_inputs)

# Spatial exploited biomass -----------------------------------------------
reg_name <- list.files(base_folder, pattern = "^(e|w)")

# Load exploited biomass estimated for all regions and resolutions
bio_data <- reg_name |> 
  map(\(x) file.path(base_folder, x, "gridded_dbpm_outputs")) |> 
  map(\(x) list.files(x, "mean_pred_det_10g-1t-bio_exploit-bio_", 
                      full.names = T, recursive = T)) |> 
  unlist() |> 
  map(\(x) read_parquet(x)) |> 
  bind_rows() |> 
  # If estimated biomass is negative, change to zero
  # If estimated biomass is larger than 1e3, change to NA
  mutate(tot_expl_pred_bio = case_when(tot_expl_pred_bio < 0 ~ 0, 
                                       tot_expl_pred_bio > 1e3 ~ NA,
                                       .default = tot_expl_pred_bio),
         tot_expl_det_bio = case_when(tot_expl_det_bio < 0 ~ 0, 
                                      tot_expl_det_bio > 1e3 ~ NA,
                                      .default = tot_expl_det_bio)) |> 
  # Add biomass for detritivores and predators
  mutate(tot_expl_bio = tot_expl_pred_bio+tot_expl_det_bio, 
         year = year(time)) |> 
  select(year, region, resolution, tot_expl_bio) |> 
  # Calculate mean total exploited biomass per year, resolution and region
  group_by(year, region, resolution) |>
  summarise(mean_expl_bio = mean(tot_expl_bio, na.rm = T)) 
  
# Merge exploited biomass for gridded and non-gridded DBPM
bio_data <- exp_bio_data |> 
  select(names(bio_data)) |> 
  bind_rows(bio_data)
  
# Plot exploited biomass
bio_data |> 
  ggplot(aes(year, mean_expl_bio, color = region, linetype = resolution))+
  geom_line()+
  theme_bw()

# Save results
bio_data |> 
  write_parquet("outputs/mean_yr_tot_exploited_biomass_1841-2010.parquet")

# Remove variables not needed
# rm(exp_bio, fishing_params, params)


# Non-spatial estimated catches -------------------------------------------
catch_dbpm_nonspatial <- reg_name |> 
  map(\(x) file.path(base_folder, x, "run_nonspatial/1deg")) |> 
  map(\(x) list.files(x, full.names = T)) |> 
  map(\(x) read_parquet(x)) |> 
  bind_rows() |> 
  group_by(year, region) |> 
  summarise(mean_catch = mean(total_catch, na.rm = T)) |> 
  mutate(resolution = "non-spatial", .after = region)


# Spatial estimated catches -----------------------------------------------
catch_files <- reg_name |> 
  map(\(x) file.path(base_folder, x, "gridded_dbpm_outputs")) |> 
  map(\(x) list.files(x, "mean_year_catch_dbpm_[0|1]",
                      full.names = T, recursive = T)) |> 
  unlist()

catch_data <- data.frame()
for(cf in catch_files){
  fao <- str_to_upper(str_replace(str_extract(cf, "(fao-[0-9]{2})"), "-", " "))
  res <- str_extract(cf, "[0-9]{0,3}deg")
  catch_df <- read_parquet(cf, col_select = !units) |> 
    mutate(mean_catch = case_when(mean_catch < 0 ~ 0, 
                                  mean_catch > 1e3 ~ NA,
                                  .default = mean_catch)) |>  
    mutate(region = fao, resolution = res, .after = year)
  catch_data <- catch_data |> 
    bind_rows(catch_df)
}

# Merge exploited biomass for gridded and non-gridded DBPM
# bio_data <- 
catch_data <- catch_dbpm_nonspatial |> 
  bind_rows(catch_data)

# Plotting data
catch_data |> 
  ggplot(aes(year, mean_catch, color = region, linetype = resolution))+
  geom_line()+
  theme_bw()

# Save results
catch_data |> 
  write_parquet("outputs/mean_yr_estimated_catches_1841-2010.parquet")



# Calculating fishing pressure plots --------------------------------------
fishing_pressure <- catch_data |> 
  left_join(bio_data, by = c("year", "region", "resolution")) |> 
  mutate(fp = mean_catch/mean_expl_bio,
         period = ifelse(year <= 1960, "steady", "historical"))

# Save results
fishing_pressure |> 
  write_parquet("outputs/mean_fishing_pressure_yr_1841-2010.parquet")
  

# Plotting fishing pressure -----------------------------------------------
fishing_pressure |> 
  filter(period == "steady") |> 
  ggplot(aes(year, fp, color = resolution))+
  geom_line(aes(size = resolution))+
  scale_size_manual(values = c(1.6, 0.8, 0.4))+
  geom_point(aes(shape = resolution))+
  scale_color_manual(values = c("#d7301f", "#fc8d59", "#fdcc8a"))+
  facet_grid(region~., scales = "free")+
  scale_y_continuous(labels = scales::label_scientific())+
  theme_bw()

ggsave("outputs/fishing_pressure_1841-1960.tif")


fishing_pressure |> 
  filter(period == "historical") |> 
  ggplot(aes(year, fp, color = resolution))+
  geom_line(aes(size = resolution))+
  scale_size_manual("DBPM resolution", values = c(1.6, 0.8, 0.4))+
  geom_point(aes(shape = resolution))+
  scale_color_manual("DBPM resolution",
                     values = c("#d7301f", "#fc8d59", "#fdcc8a"))+
  guides(shape = guide_legend(title = "DBPM resolution"))+
  facet_grid(region~., scales = "free")+
  theme_bw()+
  labs(y = "Fishing pressure")+
  theme(strip.text = element_text(family = "sans", size = 12),
        axis.text = element_text(family = "sans", size = 12), 
        axis.title.x = element_blank(), legend.position = "top",
        axis.title.y = element_text(family = "sans", size = 14, 
                                    margin = margin(5, 5, 7.5, 5, unit = "pt")), 
        legend.direction = "horizontal", legend.title.position = "top",
        legend.title = element_text(family = "sans", face = "bold", 
                                    hjust = 0.5), 
        legend.text = element_text(family = "sans", size = 12), 
        panel.grid.minor = element_blank(), 
        plot.margin = margin(0, 5, 5, 5, unit = "pt"), 
        legend.margin = margin(5, 5, 0, 5, unit = "pt"))
  
ggsave("outputs/fishing_pressure_1961-2010.png")



# Loading estimated catches forced with CCAMLR effort ---------------------
# Gridded DBPM
catch_files <- reg_name |> 
  map(\(x) file.path(base_folder, x, "gridded_dbpm_outputs")) |> 
  map(\(x) list.files(x, "mean_year_catch_dbpm_c",
                      full.names = T, recursive = T)) |> 
  unlist()

catch_ccamlr_data <- data.frame()
for(cf in catch_files){
  fao <- str_to_upper(str_replace(str_extract(cf, "(fao-[0-9]{2})"), "-", " "))
  res <- str_extract(cf, "[0-9]{0,3}deg")
  catch_df <- read_parquet(cf, col_select = !units) |> 
    mutate(mean_catch = case_when(mean_catch < 0 ~ 0, 
                                  mean_catch > 1e3 ~ NA,
                                  .default = mean_catch)) |>  
    mutate(region = fao, resolution = res, .after = year)
  catch_ccamlr_data <- catch_ccamlr_data |> 
    bind_rows(catch_df)
}

# Non-gridded DBPM
catch_dbpm_nonspatial <- reg_name |> 
  map(\(x) file.path(base_folder, x, "run_nonspatial")) |> 
  map(\(x) list.files(x, "dbpm_non", full.names = T)) |> 
  map(\(x) read_parquet(x)) |> 
  bind_rows() |> 
  filter(year >= 1969) |> 
  group_by(year, region) |> 
  summarise(mean_catch = mean(total_catch, na.rm = T)) |> 
  mutate(resolution = "non-spatial", .after = region)

# Join gridded and non-gridded data frames
catch_ccamlr_data <- catch_dbpm_nonspatial |> 
  bind_rows(catch_ccamlr_data)

catch_ccamlr_data |> 
  ggplot(aes(year, mean_catch, color = resolution))+
  geom_line()+
  facet_grid(region~., scales = "free")

# Save results
catch_ccamlr_data |> 
  write_parquet("outputs/mean_yr_estimated_catches_1969-2010.parquet")






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





base_dir <- "/g/data/vf71/la6889/dbpm_inputs"
reg_name <- list.files(base_dir, pattern = "^(e|w)")

clim_data <- reg_name |> 
  map_chr(\(x) file.path(base_dir, x, "monthly_weighted", "1deg")) |> 
  list.files(pattern = "dbpm_clim-fish-inputs_", full.names = T) |> 
  map(~read_parquet(., col_select = c(region, time, year, lphy, sphy, 
                                      tob, tos))) |> 
  bind_rows() |> 
  filter(year >= 1961) |> 
  group_by(region, year) |> 
  summarise(across(lphy:tos, 
                   #Listing statistics to be calculated
                   list(mean = ~mean(.x, na.rm = T)), 
                   #Setting column names
                   .names = "{.col}")) |> 
  ungroup() |> 
  rowwise() |> 
  mutate(phyc_vint = sphy+lphy) |> 
  select(!lphy)

# clim_data_025deg <- reg_name |>
#   map_chr(\(x) file.path(base_dir, x, "monthly_weighted", "025deg")) |>
#   list.files(pattern = "dbpm_clim-fish-inputs_", full.names = T) |>
#   map(~read_parquet(., col_select = c(region, time, year, lphy, sphy,
#                                       tob, tos))) |>
#   bind_rows() |>
#   filter(year >= 1961) |>
#   group_by(region, year) |>
#   summarise(across(lphy:tos,
#                    #Listing statistics to be calculated
#                    list(mean = ~mean(.x, na.rm = T)),
#                    #Setting column names
#                    .names = "{.col}")) |>
#   ungroup() |>
#   mutate(res = "025deg",
#          phyc_vint = sphy+lphy, .before = sphy) |>
#   rename("phypico_vint" = "sphy") |>
#   select(!lphy)
# 
# clim_data <- clim_data_1deg |>
#   bind_rows(clim_data_025deg)

# fig_temp <-
clim_data |> 
  ggplot()+
  geom_line(aes(x = year, y = tos, color = "Sea surface temperature"))+
  geom_line(aes(x = year, y = tob, 
                color = "Sea potential temperature at seafloor"))+
  facet_wrap(~region, ncol = 2, scales = "free")+
  scale_color_manual(values = c("#1f78b4", "#33a02c"),
                     labels = c("Sea potential temperature at seafloor",
                                "Sea surface temperature"))+
  labs(y = "Temperature (°C)", subtitle = "Area-weighted annual means")+
  guides(color = guide_legend(title = "DBPM inputs"))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(family = "sans", size = 14),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(family = "sans", size = 12),
        axis.text = element_text(family = "sans", size = 12),
        plot.subtitle = element_text(family = "sans", size = 12, hjust = 0.5,
                                     face = "bold"),
        legend.position = "inside", legend.position.inside = c(0.725, 0.4))

fn_out <- "composite_fig_temp_dbpm_inputs.png"
ggsave(fn_out)


clim_data |> 
  ggplot()+
  geom_line(aes(x = year, y = sphy, 
                color = "Picophytoplankton carbon concentration"))+
  geom_line(aes(x = year, y = phyc_vint,
                color = "Phytoplankton carbon concentration"))+
  facet_wrap(~region, ncol = 2, scales = "free")+
  scale_color_manual(values = c("#1f78b4", "#33a02c"),
                     labels = c("Phytoplankton carbon concentration",
                                "Picophytoplankton carbon concentration"))+
  labs(y = ~paste("Carbon concentration (mol*", m^-2, ")"),
       subtitle = "Area-weighted annual means")+
  guides(color = guide_legend(title = "DBPM inputs"))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(family = "sans", size = 14),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(family = "sans", size = 12),
        axis.text = element_text(family = "sans", size = 12),
        plot.subtitle = element_text(family = "sans", size = 12, hjust = 0.5,
                                     face = "bold"),
        legend.position = "inside", legend.position.inside = c(0.725, 0.4))

fn_out <- "composite_fig_phyto_dbpm_inputs.png"
ggsave(fn_out)



pred_det_bio <- reg_name |> 
  map_chr(\(x) file.path(base_dir, x, "gridded_dbpm_outputs")) |> 
  list.files(pattern = "10g-1t", full.names = T, recursive = T) |> 
  map(~read_parquet(.)) |> 
  bind_rows() |>
  mutate(year = year(time)) |> 
  group_by(region, resolution, year) |> 
  summarise(across(tot_expl_pred_bio:tot_det_bio, 
                   #Listing statistics to be calculated
                   list(mean = ~mean(.x, na.rm = T)), 
                   #Setting column names
                   .names = "{.col}"))

non_spat <- reg_name |> 
  map_chr(\(x) file.path(base_dir, x, "run_nonspatial")) |> 
  list.files(pattern = "10g-1t", full.names = T) |> 
  map(~read_parquet(.)) |> 
  bind_rows() |> 
  select(region, year, total_pred_biomass, total_detritivore_biomass) |> 
  group_by(region, year) |>
  summarise(across(total_pred_biomass:total_detritivore_biomass, 
                   #Listing statistics to be calculated
                   list(mean = ~mean(.x, na.rm = T)), 
                   #Setting column names
                   .names = "{.col}")) |> 
  mutate(resolution = "non-spatial", .after = region) |> 
  rename(tot_pred_bio = total_pred_biomass,
         tot_det_bio = total_detritivore_biomass)

pred_det_bio <- pred_det_bio |> 
  bind_rows(non_spat) 

pred_det_bio |> 
  ggplot()+
  geom_line(aes(x = year, y = tot_pred_bio, colour = resolution))+
  facet_wrap(~region, ncol = 2, scales = "free")+
  scale_color_manual(values = c("#d7301f", "#fc8d59", "#fdcc8a"),
                     labels = c("DBPM 0.25°", "DBPM 1°", "DBPM non-spatial"))+
  labs(y = ~paste("Biomass (g*", m^-2, ")"),
       subtitle = "Total predator biomass")+
  guides(colour = guide_legend(title = "DBPM resolution"))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(family = "sans", size = 14),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(family = "sans", size = 12),
        axis.text = element_text(family = "sans", size = 12),
        plot.subtitle = element_text(family = "sans", size = 12, hjust = 0.5,
                                     face = "bold"),
        legend.position = "inside", legend.position.inside = c(0.65, 0.375))


fn_out <- "composite_fig_total_pred_bio.png"
ggsave(fn_out)


pred_det_bio |> 
  ggplot()+
  geom_line(aes(x = year, y = tot_det_bio, colour = resolution))+
  facet_wrap(~region, ncol = 2, scales = "free")+
  scale_color_manual(values = c("#d7301f", "#fc8d59", "#fdcc8a"),
                     labels = c("DBPM 0.25°", "DBPM 1°", "DBPM non-spatial"))+
  labs(y = ~paste("Biomass (g*", m^-2, ")"),
       subtitle = "Total detritivore biomass")+
  guides(colour = guide_legend(title = "DBPM resolution"))+
  theme_bw()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(family = "sans", size = 14),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(family = "sans", size = 12),
        axis.text = element_text(family = "sans", size = 12),
        plot.subtitle = element_text(family = "sans", size = 12, hjust = 0.5,
                                     face = "bold"),
        legend.position = "inside", legend.position.inside = c(0.65, 0.375))


fn_out <- "composite_fig_total_det_bio.png"
ggsave(fn_out)



