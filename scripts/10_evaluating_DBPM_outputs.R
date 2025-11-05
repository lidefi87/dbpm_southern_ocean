# Evaluating DBPM estimated catches against observations using the observation
# range adjusted (ORA) method

# Loading libraries -------------------------------------------------------
library(dplyr)
library(ggplot2)
library(arrow)
library(tidyr)
library(cowplot)
library(stringr)
library(grid)
library(gridExtra)

# Define function to load catches data ------------------------------------
get_catch_dbpm_obs <- function(base_dir){
  catch_dbpm_nonspatial <- list.files(
    file.path(base_dir, "run_nonspatial/1deg"), full.names = T) |> 
    read_parquet() |> 
    select(c(time, year, total_catch)) |> 
    group_by(year) |> 
    summarise(catch_calib_org = mean(total_catch, na.rm = T)) #|> 
    # mutate(run = "original", res = "calibration")
  
  catch_dbpm_1deg <- list.files(
    file.path(base_dir, "gridded_dbpm_outputs/1deg"), 
    pattern = "mean_year_catch_dbpm_1", full.names = T) |> 
    read_parquet(col_select = !units) |> 
    rename(catch_1_org = mean_catch) #|>
    # mutate(run = "original", res = "1deg")
  
  catch_si_dbpm_1deg <- list.files(
    file.path(base_dir, "gridded_dbpm_outputs/1deg"),
    pattern = "mean_year_catch_dbpm_simask", full.names = T) |> 
    read_parquet(col_select = !units) |> 
    rename(catch_1_si = mean_catch) #|>
    # mutate(run = "simask", res = "1deg")
  
  catch_dbpm <- list.files(
    file.path(base_dir, "gridded_dbpm_outputs/025deg"),
    pattern = "mean_year_catch_dbpm_025", full.names = T) |> 
    read_parquet(col_select = !units) |> 
    rename(catch_025_org = mean_catch) |>
    # mutate(run = "original", res = "025deg") |> 
    # bind_rows(catch_dbpm_1deg, catch_si_dbpm_1deg, catch_dbpm_nonspatial) |> 
    left_join(catch_dbpm_1deg, by = "year") |>
    left_join(catch_si_dbpm_1deg, by = "year") |>
    left_join(catch_dbpm_nonspatial, by = "year") |>
    filter(year >= 1961) |> 
    pivot_longer(cols = starts_with("catch"), names_to = c("res", "run"), 
                 names_prefix = "catch_", names_pattern = "(.*)_(.*)", 
                 values_to = "vals")
  
  catch_obs <- list.files(file.path(base_dir, "monthly_weighted/1deg"),
                          pattern = "^dbpm_clim", full.names = T) |> 
    read_parquet() |> 
    filter(year >= 1961) 
  
  reg_code <- catch_obs |> 
    distinct(region) |> 
    pull()
  
  catch_obs <-  catch_obs|> 
    group_by(year) |> 
    summarise(obs_ccamlr = mean(catch_ccamlr*1e6, na.rm = T),
              obs_pauly = mean(catch_pauly*1e6, na.rm = T),
              obs_watson = mean(catch_tonnes_area_m2*1e6, na.rm = T)) |> 
    rowwise() |>
    mutate(min_catch_density = min(obs_ccamlr, obs_pauly, obs_watson, na.rm = T),
           max_catch_density = max(obs_ccamlr, obs_pauly, obs_watson, 
                                   na.rm = T)) |> 
    pivot_longer(cols = starts_with("obs"), names_to = "source", 
                 values_to = "obs") |> 
    mutate(min_catch_density = case_when(is.infinite(min_catch_density) ~ NA, 
                                         T ~ min_catch_density),
           max_catch_density = case_when(is.infinite(max_catch_density) ~ NA, 
                                         T ~ max_catch_density))
  
  catch_dbpm <- catch_dbpm |>
    left_join(catch_obs, by = "year") |> 
    # bind_cols(catch_obs |> 
    #             select(!year)) |> 
    #Calculating pseudo observations
    mutate(pseudo = case_when(vals < min_catch_density ~ min_catch_density,
                              vals > max_catch_density ~ max_catch_density,
                              .default = vals)) |> 
    mutate(region = reg_code, .after = run)
  
  return(catch_dbpm)

}


# Loading data ------------------------------------------------------------
base_dirs <- dir("/g/data/vf71/la6889/dbpm_inputs", pattern = "^[e|w]",
                 full.names = T)

catch_dbpm_obs <- data.frame()

for(fp in base_dirs){
  reg_name <- basename(fp) |> 
    str_replace_all("_", " ") |> 
    str_to_title()
  da <- get_catch_dbpm_obs(fp) |> 
    mutate(reg_name = reg_name, .after = region)
  catch_dbpm_obs <- catch_dbpm_obs |> 
    bind_rows(da)
}

# Save results
catch_dbpm_obs |> 
  write_parquet("outputs/simulated_observed_catches_all.parquet")
# catch_dbpm_obs <- catch_dbpm_obs |>
#   group_split(region)

# Plotting modelled vs observed catches  ----------------------------------
si_catch <- catch_dbpm_obs |> 
  #Select runs with sea ice mask and calibration (non-spatial)
  filter(run == "si" | res == "calib") |> 
  ggplot()+
  geom_ribbon(aes(x = year, ymin = min_catch_density, 
                  ymax = max_catch_density), fill = "#d8d6dd")+
  geom_line(aes(year, vals, color = res), linewidth = 0.9)+
  geom_point(aes(year, vals, color = res, shape = res), size = 1.5)+
  geom_line(aes(year, obs, color = source), linetype = "dashed")+
  geom_point(aes(year, obs, color = source, shape = source), size = 0.2)+
  facet_grid(region~., scales = "free_y")+
  scale_color_manual(values = c(#"#d7301f",
    "#fc8d59", "#fdcc8a",
    "#045a8d", "#2b8cbe", "#74a9cf"),
    labels = c(#"0.25°", 
      "1°", "calibration",
      "CCAMLR", "Pauly & Zeller", "Watson & Tidd"))+
  scale_shape_discrete(labels = c(#"0.25°", 
    "1°", "calibration",
    "CCAMLR", "Pauly & Zeller", "Watson & Tidd"))+
  # scale_linetype_manual(values = c(1, 6), 
  #                       labels = c("original", "sea ice mask"))+
  scale_x_continuous(breaks = as.integer(seq(1960, 2010, by = 10)))+
  ylab(~paste("Mean annual catch per unit area (t*", km^-2, ")"))+
  # ylab(expression(atop("Mean annual catch per unit area", 
  #                      ~paste("(t ", km^-2, ")"))))+
  # labs(subtitle = paste0(name, " (", code, ")"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(family = "sans", size = 14),
        axis.text = element_text(family = "sans", size = 12),
        legend.title = element_blank(), 
        legend.byrow = T, legend.direction = "horizontal", 
        strip.text = element_text(family = "sans", size = 12),
        legend.position = "bottom",
        legend.text = element_text(family = "sans", size = 12), 
        plot.margin = margin(5, 15, 5, 5, "pt"),
        plot.subtitle = element_text(family = "sans", size = 12, hjust = 0.5))


#Comparison original and sea ice mask runs
catch_dbpm_obs |> 
  ggplot()+
  geom_ribbon(aes(x = year, ymin = min_catch_density, 
                  ymax = max_catch_density), fill = "#d8d6dd")+
  geom_line(aes(year, vals, color = res, linetype = run), linewidth = 0.9)+
  geom_point(aes(year, vals, color = res, shape = res), size = 1.5)+
  geom_line(aes(year, obs, color = source), linetype = "dashed")+
  geom_point(aes(year, obs, color = source, shape = source), size = 0.2)+
  facet_grid(region~., scales = "free_y")+
  scale_color_manual("", values = c("#d7301f", "#fc8d59", "#fdcc8a",
                                "#045a8d", "#2b8cbe", "#74a9cf"),
                     labels = c("0.25°", "1°", "calibration", "CCAMLR", 
                                "Pauly & Zeller", "Watson & Tidd"))+
  scale_shape_discrete("", labels = c("0.25°", "1°", "calibration", "CCAMLR",
                                  "Pauly & Zeller", "Watson & Tidd"))+
  scale_linetype_manual("Runs", values = c(1, 6), 
                        labels = c("original", "sea ice mask"))+
  scale_x_continuous(breaks = as.integer(seq(1960, 2010, by = 10)))+
  ylab(~paste("Mean annual catch per unit area (t*", km^-2, ")"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(family = "sans", size = 14),
        axis.text = element_text(family = "sans", size = 12),
        legend.title = element_blank(), 
        legend.byrow = T, legend.direction = "horizontal", 
        strip.text = element_text(family = "sans", size = 12),
        legend.position = "bottom",
        legend.text = element_text(family = "sans", size = 12), 
        plot.margin = margin(5, 15, 5, 5, "pt"),
        plot.subtitle = element_text(family = "sans", size = 12, hjust = 0.5))

# Creating plots to get legeds
# Observed catches
x <- catch_dbpm_obs |>
  ggplot()+
  geom_ribbon(aes(x = year, ymin = min_catch_density, ymax = max_catch_density),
              fill = "#f1eef6")+
  geom_line(aes(year, obs, color = source), linetype = "dashed")+
  geom_point(aes(year, obs, color = source, shape = source), size = 0.2)+
  scale_color_manual(values = c("#045a8d", "#2b8cbe", "#74a9cf"),
                     labels = c("CCAMLR", "Pauly & Zeller", "Watson & Tidd"))+
  scale_shape_discrete(labels = c("CCAMLR", "Pauly & Zeller", "Watson & Tidd"))+
  guides(color = guide_legend(title = "Observed catches"))+
  guides(shape = guide_legend(title = "Observed catches"))+
  theme_bw()+
  theme(legend.title = element_text(face = "bold"), 
        legend.byrow = T, legend.position = "right",
        legend.text = element_text(family = "sans", size = 12))

leg_obs <- get_legend(x)

# Simulated catches
y <- catch_dbpm_obs |>
  ggplot()+
  geom_ribbon(aes(x = year, ymin = min_catch_density, ymax = max_catch_density),
              fill = "#f1eef6")+
  geom_line(aes(year, vals, color = res), linewidth = 0.9)+
  geom_point(aes(year, vals, color = res, shape = res), size = 1.5)+
  scale_color_manual(values = c("#d7301f", "#fc8d59", "#fdcc8a"),
                     labels = c("0.25°", "1°", "calibration (non-spatial)"))+
  scale_shape_discrete(labels = c("0.25°", "1°", "calibration (non-spatial)"))+
  theme_bw()+
  guides(color = guide_legend(title = "DBPM estimated catches"))+
  guides(shape = guide_legend(title = "DBPM estimated catches"))+
  theme(legend.title = element_text(face = "bold"), 
        legend.byrow = T, legend.position = "right",
        legend.text = element_text(family = "sans", size = 12)) 
        
leg_est <- get_legend(y)

si_catch_noleg <- si_catch+theme(legend.position = "none")

fig <- grid.arrange(arrangeGrob(
  plot_grid(si_catch_noleg, 
            plot_grid(leg_est, leg_obs, ncol = 2, align = "hv"), 
            nrow = 2, rel_heights = c(0.8, 0.2))))

ggsave("outputs/composite_fig_catches_obs.png", fig, height = 6, width = 8)



# Calculating ORA metrics -------------------------------------------------
ora_metrics <- catch_dbpm_obs |> 
  group_by(res, run, region) |>
  summarise(pear_cor = cor(pseudo, vals),
            mae = sum(abs(vals-pseudo))/n(),
            ri = exp(sqrt(sum(log(pseudo/vals)^2)/n())),
            mef = ((sum((pseudo-mean(pseudo))^2)-sum((vals-mean(vals))^2))/
                     sum((pseudo-mean(pseudo))^2)))

fout <- file.path(base_dir, "gridded_dbpm_outputs", 
                  paste0("ora_metric_catches_all-res_", reg, 
                         "_1961-2010.parquet"))
ora_metrics |> 
  write_parquet(fout)
