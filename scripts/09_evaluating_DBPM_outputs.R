# Evaluating DBPM estimated catches against observations using the observation
# range adjusted (ORA) method

# Loading libraries -------------------------------------------------------
library(dplyr)
library(ggplot2)
library(arrow)
library(tidyr)
library(cowplot)
library(stringr)


# Define function to load catches data ------------------------------------
get_catch_dbpm_obs <- function(base_dir){
  catch_dbpm_nonspatial <- list.files(file.path(base_dir, "run_nonspatial/1deg"), 
                                      full.names = T) |> 
    read_parquet() |> 
    select(c(time, year, total_catch)) |> 
    group_by(year) |> 
    summarise(catch_nonspat = mean(total_catch, na.rm = T))
  
  catch_dbpm_1deg <- list.files(file.path(base_dir, "gridded_dbpm_outputs/1deg"),
                                pattern = "mean_year_catch", full.names = T) |> 
    read_parquet(col_select = !units) |> 
    rename(catch_1 = mean_catch)
  
  catch_dbpm <- list.files(file.path(base_dir, 
                                     "gridded_dbpm_outputs/025deg"),
                           pattern = "mean_year_catch", full.names = T) |> 
    read_parquet() |> 
    rename(catch_025 = mean_catch) |> 
    left_join(catch_dbpm_1deg, by = "year") |> 
    left_join(catch_dbpm_nonspatial, by = "year") |> 
    filter(year >= 1961) |> 
    pivot_longer(cols = starts_with("catch"), names_to = "res", 
                 values_to = "vals")
  
  catch_obs <- list.files(file.path(base_dir, "monthly_weighted/1deg"),
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
    pivot_longer(cols = starts_with("obs"), names_to = "source", 
                 values_to = "obs") |> 
    mutate(min_catch_density = case_when(is.infinite(min_catch_density) ~ NA, 
                                         T ~ min_catch_density),
           max_catch_density = case_when(is.infinite(max_catch_density) ~ NA, 
                                         T ~ max_catch_density))
  
  catch_dbpm <- catch_dbpm |> 
    bind_cols(catch_obs |> 
                select(!year)) |> 
    #Calculating pseudo observations
    mutate(pseudo = case_when(vals < min_catch_density ~ min_catch_density,
                              vals > max_catch_density ~ max_catch_density,
                              .default = vals))
  
  return(catch_dbpm)

}



# Loading data ------------------------------------------------------------
name_reg <- "weddell"
reg <- "fao-48"
base_dir <- file.path("/g/data/vf71/la6889/dbpm_inputs", name_reg)


# Plotting modelled vs observed catches  ----------------------------------
catch_dbpm_obs <- get_catch_dbpm_obs(base_dir)

catch_dbpm_obs |>
  ggplot()+
  geom_ribbon(aes(x = year, ymin = min_catch_density, ymax = max_catch_density),
              fill = "#f1eef6")+
  geom_line(aes(year, vals, color = res), linewidth = 0.9)+
  geom_point(aes(year, vals, color = res, shape = res), size = 2)+
  geom_line(aes(year, obs, color = source), linetype = "dashed")+
  geom_point(aes(year, obs, color = source, shape = source), size = 0.1)+
  scale_color_manual(values = c("#d7301f", "#fc8d59", "#fdcc8a", 
                                "#045a8d", "#2b8cbe", "#74a9cf"),
                     labels = c("DBPM 0.25°", "DBPM 1°", "DBPM non-spatial",
                                  "CCAMLR", "Pauly & Zeller", "Watson et al"))+
  scale_shape_discrete(labels = c("DBPM 0.25°", "DBPM 1°", "DBPM non-spatial",
                                  "CCAMLR", "Pauly & Zeller", "Watson et al"))+
  scale_x_continuous(breaks = as.integer(seq(1960, 2010, by = 10)))+
  ylab(~paste("Mean annual catch per unit area (g*", m^-2, ")"))+
  # ylab(expression(atop("Mean annual catch per unit area", 
  #                      ~paste("(g*", m^-2, ")"))))+
  labs(subtitle = paste0(str_to_title(str_replace_all(name_reg, "_", " ")),
                         " (", str_replace_all(str_to_upper(reg), "-", " "), 
                         ")"))+
  theme_bw()+
  theme(panel.grid.minor = element_blank(), axis.title.x = element_blank(), 
        axis.title.y = element_text(family = "sans", size = 14),
        axis.text = element_text(family = "sans", size = 12),
        # legend.position = "inside", legend.position.inside = c(0.15, 0.825), 
        legend.title = element_blank(), 
        legend.byrow = T, legend.direction = "horizontal", 
        legend.position = "bottom",
        legend.text = element_text(family = "sans", size = 12), 
        plot.margin = margin(c(5, 15, 5, 5), "cm"),
        plot.subtitle = element_text(family = "sans", size = 12, hjust = 0.5))

fout <- file.path(base_dir, "gridded_dbpm_outputs", 
                  paste0("mean_catches_ts_all-res_", reg, "_1961-2010.png"))

ggsave(fout, device = "png")



# Calculating ORA metrics -------------------------------------------------
ora_metrics <- catch_dbpm_obs |> 
  group_by(res) |>
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
