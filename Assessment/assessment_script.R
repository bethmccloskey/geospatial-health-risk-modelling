# ==============================================================================
# Spatial Epidemiology in Public Health: Côte d'Ivoire (2021 DHS)
# Script: assessment_script.R
# Author: Beth (MSc Spatial Epidemiology)
# 
# Description:
#   Full end-to-end analysis script covering:
#   1. Exploratory Spatial Data Analysis (ESDA) of Open Defecation (Admin 2)
#      - Spatial weights construction (Queen contiguity, row-standardized)
#      - Global Moran's I & 9,999 Monte Carlo permutations
#      - Local Moran's I (LISA) hotspot and outlier classification
#   2. Geostatistical Modeling of Childhood Malaria (RDT cluster points)
#      - Empirical logit variance-stabilizing transformation
#      - Empirical semivariogram with 95% Monte Carlo CSR envelope
#      - Theoretical semivariogram modeling (Exponential fit via WLS)
#      - Ordinary Kriging spatial interpolation on a 5 km national grid
#      - 10-fold internal cross-validation
#   3. Population-at-Risk Analysis (WorldPop 1 km gridded population)
#      - High-risk masking (80th percentile prevalence threshold)
#      - National and district-level zonal population estimation
#      - Kriging prediction uncertainty sensitivity analysis
#      - Exporting spatial layers (GeoTIFF, GeoPackage, GeoJSON) for QGIS
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Setup & Dependencies
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(janitor)
  library(scales)
  library(ggplot2)
  library(gt)
  library(tableone)
  
  # Spatial vector and raster processing
  library(sf)
  library(terra)
  library(tmap)
  
  # Spatial econometrics and geostatistics
  library(spdep)
  library(sp)
  library(gstat)
  library(automap)
  library(geoR)
})

# Source custom geostatistical helper functions
if (file.exists("functions.R")) {
  source("functions.R")
} else {
  stop("Helper functions file 'functions.R' not found.")
}

# ------------------------------------------------------------------------------
# 2. Data Ingestion & Validation
# ------------------------------------------------------------------------------
cat("Loading survey datasets and spatial boundaries...\n")

# A. Malaria RDT survey clusters (n = 468)
adm1_rdt <- read.csv("data/cluster_rdt.csv")

# B. Open defecation aggregated to Admin 2 departments (n = 106)
adm2_def <- read.csv("data/adm2_opendef.csv")

# C. Spatial boundaries: Admin 2 DHS shapefile
adm2_dhs_def <- st_read("data/civ_DHS.shp", quiet = TRUE) %>%
  clean_names()

# D. National (Admin 0) and Regional/District (Admin 1) boundaries
adm0 <- st_read("data/civ_adm0.shp", quiet = TRUE)
adm1 <- st_read("data/civ_adm1.shp", quiet = TRUE)

# E. 2021 WorldPop 1 km constrained population raster
civ_pop <- rast("data/civ_pop_2021_CN_1km_R2025A_UA_v1.tif")

cat("Data successfully loaded.\n")
cat(" - Malaria RDT Clusters:", nrow(adm1_rdt), "\n")
cat(" - Open Defecation Departments:", nrow(adm2_def), "\n")

# ==============================================================================
# SECTION 1: Exploratory Spatial Data Analysis of Open Defecation
# ==============================================================================

cat("\n--- Running Section 1: Open Defecation Spatial Analysis ---\n")

# 1.1 Data Harmonization
# Join tabular survey responses to department boundaries
full_def <- adm2_dhs_def %>%
  inner_join(adm2_def, by = c("adm2_pcode", "dhs_region")) %>%
  select(-n_clusters, -n_househol, -n_opendef)

# Compute district-level urban proportion from cluster data
urban_prop <- adm1_rdt %>%
  group_by(adm1name) %>%
  summarise(
    urban_prop = mean(urban, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(adm1name = tolower(adm1name))

# Append urban proportion to departmental spatial data
full_def <- full_def %>%
  left_join(urban_prop, by = c("dhs_region" = "adm1name"))

# 1.2 Summary Statistics: District Table
regional_summary <- full_def %>%
  st_drop_geometry() %>%
  group_by(dhs_region) %>%
  summarise(
    no_households = sum(n_hhold),
    no_open_def   = sum(n_open),
    avg_prop_def  = round(mean(prop_od) * 100, 2),
    urban_prop    = round(mean(urban_prop) * 100, 2),
    .groups       = "drop"
  ) %>%
  mutate(dhs_region = str_to_title(dhs_region))

cat("Regional summary calculated across", nrow(regional_summary), "districts.\n")

# 1.3 Exploratory Visualizations
# A. Histogram of department open defecation proportions
dist_alone <- ggplot(full_def, aes(x = prop_od)) +
  geom_histogram(binwidth = 0.05, fill = "#8B1E3F", colour = "white", alpha = 0.9) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  labs(
    title = "Distribution of Department Open Defecation Proportions",
    subtitle = "Right-skewed distribution across 106 departments",
    x = "Proportion of Households Practicing Open Defecation",
    y = "Frequency"
  ) +
  theme_classic(base_size = 14)

# B. Scatter plot: Urbanization vs. Open Defecation
region_scat <- ggplot(regional_summary, aes(x = urban_prop, y = avg_prop_def)) +
  geom_point(size = 3.5, color = "#2B5B84", alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, color = "#C0392B", linewidth = 1.2) +
  scale_x_continuous(limits = c(0, 100), labels = function(x) paste0(x, "%")) +
  scale_y_continuous(limits = c(0, 60), labels = function(x) paste0(x, "%")) +
  labs(
    title = "Urban Coverage vs. Open Defecation Prevalence",
    subtitle = "Strong negative ecological association across districts",
    x = "Average Urban Coverage (%)",
    y = "Average Open Defecation Prevalence (%)"
  ) +
  theme_classic(base_size = 14)

# 1.4 Spatial Weights Matrix Construction
# Queen contiguity: departments sharing a common border or point vertex
neighbours <- poly2nb(full_def, queen = TRUE)

# Row-standardized weights matrix (style = 'W')
weights_queen <- nb2listw(neighbours, style = "W", zero.policy = TRUE)

# 1.5 Global Spatial Autocorrelation (Global Moran's I)
def_moran <- moran.test(full_def$prop_od, listw = weights_queen, zero.policy = TRUE)
cat("Global Moran's I Statistic:", round(def_moran$estimate[1], 4), 
    "| p-value:", format.pval(def_moran$p.value), "\n")

# Standardized z-score and spatial lag
full_def <- full_def %>%
  mutate(
    scale_prop_od = as.numeric(scale(prop_od)),
    lag_scaled_prop_od = lag.listw(weights_queen, scale_prop_od, zero.policy = TRUE)
  )

# Monte Carlo permutation test (9,999 simulations under CSR)
set.seed(42)
mc_moran <- moran.mc(full_def$prop_od, listw = weights_queen, nsim = 9999, zero.policy = TRUE)
cat("Monte Carlo Pseudo p-value:", mc_moran$p.value, "\n")

# Global Moran's I scatterplot
stat_label <- paste0("Global Moran's I = ", round(def_moran$estimate[1], 3), "\n(p < 0.001)")
global_moran_scat <- ggplot(full_def, aes(x = scale_prop_od, y = lag_scaled_prop_od)) +
  geom_point(color = "#1F4E79", size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", color = "#C0392B", se = FALSE, linewidth = 1.1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  annotate("text", x = 1.2, y = -1.2, label = stat_label, hjust = 0, size = 5, fontface = "bold", color = "#C0392B") +
  labs(
    title = "Moran's I Scatterplot of Open Defecation in Côte d'Ivoire",
    subtitle = "Positive spatial autocorrelation indicates strong geographical clustering",
    x = "Standardised Open Defecation Proportion (Z-score)",
    y = "Spatially Lagged Average Z-score (Wz)"
  ) +
  theme_minimal(base_size = 14)

# 1.6 Local Indicators of Spatial Association (LISA)
local_moran <- localmoran(full_def$prop_od, weights_queen, zero.policy = TRUE, alternative = "two.sided")
colnames(local_moran) <- c("local_I", "expected_I", "variance_I", "z_statistic", "p_value")
full_def$p_value <- local_moran[, "p_value"]

# LISA classification function
classify_lisa <- function(data, alpha = 0.05) {
  data %>%
    mutate(
      cluster = case_when(
        p_value > alpha ~ "Not Significant",
        scale_prop_od > 0 & lag_scaled_prop_od > 0 ~ "High-High",
        scale_prop_od <= 0 & lag_scaled_prop_od <= 0 ~ "Low-Low",
        scale_prop_od > 0 & lag_scaled_prop_od <= 0 ~ "High-Low",
        scale_prop_od <= 0 & lag_scaled_prop_od > 0 ~ "Low-High"
      ),
      cluster = factor(cluster, levels = c("High-High", "Low-Low", "High-Low", "Low-High", "Not Significant"))
    )
}

lisa_05 <- classify_lisa(full_def, alpha = 0.05)
lisa_01 <- classify_lisa(full_def, alpha = 0.01)

# Export spatial layers for QGIS cartography
st_write(lisa_05, "morans_i_local2_05.geojson", delete_dsn = TRUE, quiet = TRUE)
st_write(lisa_01, "morans_i_local2_01.geojson", delete_dsn = TRUE, quiet = TRUE)
cat("Section 1 spatial vectors successfully exported to GeoJSON.\n")

# ==============================================================================
# SECTION 2: Geostatistical Modeling of Childhood Malaria
# ==============================================================================

cat("\n--- Running Section 2: Malaria Geostatistics & Kriging ---\n")

# 2.1 Calculate Cluster Prevalence & Empirical Logit Transformation
adm1_rdt <- adm1_rdt %>%
  mutate(
    rdt_prev_child = round(n_rdt_pos / n_child, 3),
    # Empirical logit transformation with +0.5 continuity correction
    logitp = log((n_rdt_pos + 0.5) / (n_child - n_rdt_pos + 0.5))
  )

# Reproject cluster coordinates to EPSG:3857 (meters)
adm1_rdt_sp <- adm1_rdt %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326) %>%
  st_transform(crs = 3857)

adm1_spatial <- as(adm1_rdt_sp, "Spatial")

# 2.2 Empirical Semivariogram & Theoretical Model Fitting
cat("Fitting theoretical semivariogram via automap...\n")
variofit <- autofitVariogram(
  formula    = logitp ~ 1,
  input_data = adm1_spatial,
  model      = c("Exp", "Sph", "Gau")
)

vgm_model <- variofit$var_model
nugget_val <- vgm_model$psill[vgm_model$model == "Nug"]
partial_sill <- vgm_model$psill[vgm_model$model != "Nug"]
total_sill <- nugget_val + partial_sill
range_km <- vgm_model$range[vgm_model$model != "Nug"] / 1000
practical_range_km <- 3 * range_km

cat("Fitted Model:", as.character(vgm_model$model[vgm_model$model != "Nug"]), "\n")
cat(" - Nugget:", round(nugget_val, 3), "\n")
cat(" - Sill:", round(total_sill, 3), "\n")
cat(" - Range Parameter:", round(range_km, 2), "km\n")
cat(" - Practical Range (~3a):", round(practical_range_km, 2), "km\n")

# 2.3 Ordinary Kriging Interpolation
cat("Generating 5 km prediction grid and running Ordinary Kriging...\n")
adm0_proj <- st_transform(adm0, 3857)
study_area <- st_make_valid(adm0_proj) %>% st_union()

# Regular 5 km prediction grid
pred_grid <- st_make_grid(study_area, cellsize = 5000, what = "centers") %>%
  st_as_sf()
pred_locations <- pred_grid[st_within(pred_grid, study_area, sparse = FALSE)[, 1], ]

# Run Ordinary Kriging
orkriging <- autoKrige(
  formula    = logitp ~ 1,
  input_data = adm1_spatial,
  new_data   = as(pred_locations, "Spatial"),
  model      = c("Exp", "Sph", "Gau"),
  verbose    = FALSE
)

# Extract predictions and back-transform from logit to prevalence scale
pred_logit <- orkriging$krige_output$var1.pred
pred_se    <- orkriging$krige_output$var1.stdev

pred_prevalence <- plogis(pred_logit)
lower_ci        <- plogis(pred_logit - 1.96 * pred_se)
upper_ci        <- plogis(pred_logit + 1.96 * pred_se)
ci_width        <- upper_ci - lower_ci

coords_grid <- st_coordinates(pred_locations)
pred_df <- data.frame(
  x          = coords_grid[, 1],
  y          = coords_grid[, 2],
  prevalence = pred_prevalence,
  lower      = lower_ci,
  upper      = upper_ci,
  ci_width   = ci_width
)

# Rasterize continuous prediction surface
pred_rast <- terra::rast(pred_df, type = "xyz", crs = st_crs(adm1_rdt_sp)$wkt)
study_area_vect <- terra::vect(study_area)
pred_rast <- terra::crop(pred_rast, study_area_vect)
pred_rast <- terra::mask(pred_rast, study_area_vect)

writeRaster(pred_rast, "kriging_prevalence_ci2.tif", overwrite = TRUE)
cat("Kriging raster exported to 'kriging_prevalence_ci2.tif'.\n")

# 2.4 10-Fold Cross-Validation
cat("Evaluating predictive performance with 10-fold cross-validation...\n")
set.seed(123)
cv_res <- autoKrige.cv(
  formula    = logitp ~ 1,
  input_data = adm1_spatial,
  nfold      = 10,
  verbose    = c(FALSE, FALSE),
  model      = c("Exp", "Sph", "Gau")
)

comp_cv <- compare.cv(cv_res, col.names = c("Value"))
print(comp_cv)

# ==============================================================================
# SECTION 3: Population at Risk Using Gridded Population Data
# ==============================================================================

cat("\n--- Running Section 3: Population-at-Risk Analysis ---\n")

prev_rast  <- pred_rast[["prevalence"]]
lower_rast <- pred_rast[["lower"]]
upper_rast <- pred_rast[["upper"]]

# Project population raster to match kriged surface CRS
pop_proj <- project(civ_pop, crs(prev_rast))
prev_resamp  <- resample(prev_rast, pop_proj, method = "bilinear")
lower_resamp <- resample(lower_rast, pop_proj, method = "bilinear")
upper_resamp <- resample(upper_rast, pop_proj, method = "bilinear")

# Define 80th percentile high-risk threshold
thr_80 <- global(prev_resamp, quantile, probs = 0.8, na.rm = TRUE)[1, 1]
cat("80th percentile prevalence threshold:", round(thr_80, 4), "\n")

# High-risk mask and population at risk
highrisk_mask <- prev_resamp >= thr_80
pop_highrisk  <- pop_proj * highrisk_mask
names(pop_highrisk) <- "pop_highrisk"

nat_high_pop  <- global(pop_highrisk, "sum", na.rm = TRUE)[1, 1]
nat_total_pop <- global(pop_proj, "sum", na.rm = TRUE)[1, 1]
nat_pct_high  <- 100 * (nat_high_pop / nat_total_pop)

cat(sprintf("Total National Population: %.2f million\n", nat_total_pop / 1e6))
cat(sprintf("Population in High-Risk Zones: %.2f million (%.2f%%)\n", 
            nat_high_pop / 1e6, nat_pct_high))

# Dissolve departments into districts
adm_district <- adm2_dhs_def %>%
  group_by(dhs_region) %>%
  summarise(geometry = st_union(geometry), .groups = "drop")

adm_district_v <- project(vect(adm_district), crs(pop_proj))

adm_district$pop_total    <- extract(pop_proj, adm_district_v, fun = sum, na.rm = TRUE)[, 2]
adm_district$pop_highrisk <- extract(pop_highrisk, adm_district_v, fun = sum, na.rm = TRUE)[, 2]

adm1_tbl <- adm_district %>%
  st_drop_geometry() %>%
  mutate(
    pct_highrisk = round(100 * pop_highrisk / pop_total, 2),
    dhs_region   = str_to_title(dhs_region)
  ) %>%
  arrange(desc(pop_highrisk))

print(adm1_tbl[, c("dhs_region", "pop_highrisk", "pop_total", "pct_highrisk")])

# Uncertainty sensitivity analysis
pop_high_lower <- pop_proj * (lower_resamp >= thr_80)
pop_high_upper <- pop_proj * (upper_resamp >= thr_80)

nat_low_pop <- global(pop_high_lower, "sum", na.rm = TRUE)[1, 1]
nat_up_pop  <- global(pop_high_upper, "sum", na.rm = TRUE)[1, 1]

cat("\nSensitivity Analysis (95% CI uncertainty):\n")
cat(" - Conservative (Lower Bound >= 0.54):", round(nat_low_pop / 1e6, 2), "million\n")
cat(" - Expected (Point Prediction >= 0.54):", round(nat_high_pop / 1e6, 2), "million\n")
cat(" - Inclusive (Upper Bound >= 0.54):", round(nat_up_pop / 1e6, 2), "million\n")

# 3.4 Export all rasters and layers for QGIS cartography
cat("\nExporting aligned GeoTIFFs and GeoPackage for QGIS cartography...\n")
writeRaster(prev_resamp, "prevalence_aligned.tif", overwrite = TRUE)
writeRaster(highrisk_mask, "highrisk_mask.tif", overwrite = TRUE)
writeRaster(pop_highrisk, "pop_highrisk.tif", overwrite = TRUE)
writeRaster(pop_high_lower, "pop_highrisk_lower.tif", overwrite = TRUE)
writeRaster(pop_high_upper, "pop_highrisk_upper.tif", overwrite = TRUE)
st_write(adm_district, "adm_district.gpkg", layer = "adm_district", delete_layer = TRUE, quiet = TRUE)

cat("\nAnalysis complete! All outputs generated successfully.\n")
