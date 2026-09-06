# ============================================================
# SESGO 1
# Sampling Bias in Botanical Herbarium Collections
# Southern Ecuador
# ============================================================
#
# Version: 0.1
# Status: Under active development
#
# This script evaluates spatial and environmental sampling bias
# in botanical herbarium collections from southern Ecuador.
#
# Main analyses include data quality assessment, spatial filtering,
# accessibility bias (distance to roads), elevational bias,
# comparison with random reference points, and binomial
# generalized linear modelling.
#
# NOTE:
# This script represents the current analytical workflow and may
# be refined as the project develops.
#
# DATA AVAILABILITY:
# Raw occurrence data and some spatial datasets are not included
# in this repository. See README.md for data sources, availability,
# and download instructions.
# ============================================================

# ============================================================
# 1. LOAD LIBRARIES
# ============================================================

library(sf)
library(terra)
library(tidyverse)
library(janitor)
library(geodata)


# ============================================================
# 2. LOAD DATA
# ============================================================

datosA <- read_tsv(
  "data/raw/occurrence_A.txt",
  comment = "#",
  na = c("", "NA")
) %>%
  clean_names()

cat("Data loaded. Rows:", nrow(datosA), "\n")


# ============================================================
# 3. INITIAL DATA EXPLORATION
# ============================================================

# Dataset dimensions
dim(datosA)

# Dataset structure
glimpse(datosA)


# ------------------------------------------------------------
# Missing values by variable
# ------------------------------------------------------------

na_summary <- datosA %>%
  summarise(
    across(
      everything(),
      ~ mean(is.na(.))
    )
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "variable",
    values_to = "na_ratio"
  ) %>%
  arrange(desc(na_ratio))


# Variables with the highest proportion of missing values
head(na_summary, 20)

# Variables with the lowest proportion of missing values
tail(na_summary, 20)


# ============================================================
# 4. DATA QUALITY CHECKS
# ============================================================

# ------------------------------------------------------------
# Coordinate ranges
# ------------------------------------------------------------

summary(datosA$decimal_latitude)
summary(datosA$decimal_longitude)


# ------------------------------------------------------------
# Temporal coverage
# ------------------------------------------------------------

summary(datosA$year)

datosA %>%
  count(year) %>%
  arrange(year)


# ============================================================
# 5. DATA SELECTION AND CLEANING
# ============================================================

datosA_clean <- datosA %>%
  select(
    occurrence_id,
    decimal_latitude,
    decimal_longitude,
    year,
    basis_of_record
  ) %>%
  filter(
    !is.na(decimal_latitude),
    !is.na(decimal_longitude),

    # Basic geographic coordinate validation
    decimal_latitude >= -90,
    decimal_latitude <= 90,
    decimal_longitude >= -180,
    decimal_longitude <= 180
  ) %>%
  distinct(
    decimal_longitude,
    decimal_latitude,
    .keep_all = TRUE
  )

cat(
  "Records after coordinate cleaning:",
  nrow(datosA_clean),
  "\n"
)


if (nrow(datosA_clean) == 0) {
  stop("ERROR: No records remain after cleaning.")
}


# ============================================================
# 6. STUDY AREA
# ============================================================

cat("Loading administrative boundaries...\n")


# Administrative boundaries from GADM

Ecuador <- gadm(
  country = "ECU",
  level = 1,
  path = tempdir()
)

Ecuador_sf <- st_as_sf(Ecuador)


# Study region

region <- Ecuador_sf %>%
  filter(
    NAME_1 %in% c(
      "Azuay",
      "Cañar",
      "Loja",
      "Zamora Chinchipe",
      "Morona Santiago"
    )
  ) %>%
  st_union()


if (is.null(region) || length(region) == 0) {
  stop("ERROR: Study region could not be created.")
}


plot(
  region,
  main = "Study Region"
)


# ============================================================
# 7. CONVERT RECORDS TO SPATIAL DATA
# ============================================================

datosA_sf <- st_as_sf(
  datosA_clean,
  coords = c(
    "decimal_longitude",
    "decimal_latitude"
  ),
  crs = 4326
)


# ============================================================
# 8. SPATIAL FILTERING
# ============================================================

# Match the CRS of the occurrence records to the study region

datosA_sf <- st_transform(
  datosA_sf,
  st_crs(region)
)


# Retain records within the study region

datosA_sf <- st_intersection(
  datosA_sf,
  region
)


cat(
  "Records within study region:",
  nrow(datosA_sf),
  "\n"
)


# ------------------------------------------------------------
# Spatial visualization
# ------------------------------------------------------------

plot(
  st_geometry(region),
  col = "lightgrey",
  main = "Herbarium Records in Study Region"
)

plot(
  st_geometry(datosA_sf),
  add = TRUE,
  col = "red",
  pch = 16,
  cex = 0.6
)


# ============================================================
# 9. PROJECT TO UTM
# ============================================================

# UTM Zone 17S
# Distances and spatial measurements are calculated in metres.

crs_utm <- 32717


region <- st_transform(
  region,
  crs_utm
)

datosA_sf <- st_transform(
  datosA_sf,
  crs_utm
)


cat(
  "Spatial data projected to UTM Zone 17S (EPSG:",
  crs_utm,
  ")\n"
)


# ============================================================
# 10. ACCESSIBILITY BIAS
# DISTANCE TO ROADS
# ============================================================

cat("Loading road network...\n")


roads <- st_read(
  "data/raw/vias.shp",
  quiet = TRUE
)


if (nrow(roads) == 0) {
  stop("ERROR: Roads file is empty.")
}


# Transform roads to the same projected CRS

roads <- st_transform(
  roads,
  crs_utm
)


# Confirm that spatial data share the same CRS

stopifnot(
  st_crs(roads) == st_crs(datosA_sf)
)


# ------------------------------------------------------------
# Distance from herbarium records to roads
# ------------------------------------------------------------

cat(
  "Calculating distances from herbarium records to roads...\n"
)


dist_matrix_obs <- st_distance(
  datosA_sf,
  roads
)


datosA_sf$dist_road <- apply(
  as.matrix(dist_matrix_obs),
  1,
  min
)


summary(datosA_sf$dist_road)


# ============================================================
# 11. RANDOM REFERENCE POINTS
# ============================================================

# Set seed for reproducibility

set.seed(123)


cat(
  "Generating random reference points...\n"
)


random_points <- st_sample(
  region,
  size = 5000
)


random_points_sf <- st_as_sf(
  random_points
)


# ------------------------------------------------------------
# Distance from random points to roads
# ------------------------------------------------------------

cat(
  "Calculating distances from random points to roads...\n"
)


stopifnot(
  st_crs(random_points_sf) == st_crs(roads)
)


dist_matrix_rand <- st_distance(
  random_points_sf,
  roads
)


random_points_sf$dist_road <- apply(
  as.matrix(dist_matrix_rand),
  1,
  min
)


summary(random_points_sf$dist_road)


# ============================================================
# 12. ACCESSIBILITY COMPARISON
# ============================================================

# ------------------------------------------------------------
# Histograms
# ------------------------------------------------------------

par(
  mfrow = c(1, 2)
)


hist(
  datosA_sf$dist_road,
  breaks = 50,
  main = "Distance to Roads:\nHerbarium Records",
  xlab = "Distance (m)"
)


hist(
  random_points_sf$dist_road,
  breaks = 50,
  main = "Distance to Roads:\nRandom Points",
  xlab = "Distance (m)"
)


par(
  mfrow = c(1, 1)
)


# ------------------------------------------------------------
# Boxplot
# ------------------------------------------------------------

boxplot(
  datosA_sf$dist_road,
  random_points_sf$dist_road,
  names = c(
    "Collections",
    "Random"
  ),
  ylab = "Distance to roads (m)",
  main = "Distance to Roads Comparison"
)


# ------------------------------------------------------------
# Wilcoxon rank-sum test
# ------------------------------------------------------------

wilcox_test_roads <- wilcox.test(
  datosA_sf$dist_road,
  random_points_sf$dist_road
)


print(
  wilcox_test_roads
)


# ============================================================
# 13. DIGITAL ELEVATION MODEL
# ============================================================

cat(
  "Loading and processing DEM data...\n"
)


# ------------------------------------------------------------
# Digital elevation model tiles
# ------------------------------------------------------------

files <- list.files(
  "data/raw/DEM",
  pattern = "\\.tif$",
  full.names = TRUE
)


if (length(files) == 0) {
  stop(
    "ERROR: No DEM files found in data/raw/DEM."
  )
}


dem_list <- lapply(
  files,
  terra::rast
)


# ------------------------------------------------------------
# Mosaic DEM tiles
# ------------------------------------------------------------

dem <- do.call(
  terra::mosaic,
  dem_list
)


cat(
  "DEM mosaic created.\n"
)


# ------------------------------------------------------------
# Project DEM to UTM Zone 17S
# ------------------------------------------------------------

dem_utm <- terra::project(
  dem,
  "EPSG:32717"
)


# ------------------------------------------------------------
# Convert study region to SpatVector
# ------------------------------------------------------------

region_vect <- terra::vect(
  region
)


# ------------------------------------------------------------
# Crop and mask DEM
# ------------------------------------------------------------

dem_crop <- terra::crop(
  dem_utm,
  region_vect
)


dem_mask <- terra::mask(
  dem_crop,
  region_vect
)


cat(
  "DEM cropped and masked to study region.\n"
)


# ============================================================
# 14. EXTRACT ELEVATION
# ============================================================

# ------------------------------------------------------------
# Herbarium records
# ------------------------------------------------------------

cat(
  "Extracting elevation for herbarium records...\n"
)


elev_extract_obs <- terra::extract(
  dem_mask,
  terra::vect(datosA_sf)
)


datosA_sf$elevation <- elev_extract_obs[[2]]


# ------------------------------------------------------------
# Random points
# ------------------------------------------------------------

cat(
  "Extracting elevation for random points...\n"
)


elev_extract_rand <- terra::extract(
  dem_mask,
  terra::vect(random_points_sf)
)


random_points_sf$elevation <- elev_extract_rand[[2]]


# ------------------------------------------------------------
# Missing elevation values
# ------------------------------------------------------------

na_obs <- sum(
  is.na(datosA_sf$elevation)
)

na_rand <- sum(
  is.na(random_points_sf$elevation)
)


cat(
  "NA elevation values - Records:",
  na_obs,
  "\n"
)

cat(
  "NA elevation values - Random points:",
  na_rand,
  "\n"
)


# ------------------------------------------------------------
# Elevation summaries
# ------------------------------------------------------------

cat(
  "\nElevation summary - Herbarium records:\n"
)

summary(
  datosA_sf$elevation
)


cat(
  "\nElevation summary - Random points:\n"
)

summary(
  random_points_sf$elevation
)


# ============================================================
# 15. ELEVATIONAL BIAS
# ============================================================

# ------------------------------------------------------------
# Boxplot
# ------------------------------------------------------------

boxplot(
  datosA_sf$elevation,
  random_points_sf$elevation,
  names = c(
    "Collections",
    "Random"
  ),
  ylab = "Elevation (m)",
  main = "Elevation Comparison"
)


# ------------------------------------------------------------
# Wilcoxon rank-sum test
# ------------------------------------------------------------

wilcox_test_elev <- wilcox.test(
  datosA_sf$elevation,
  random_points_sf$elevation
)


print(
  wilcox_test_elev
)


# ============================================================
# 16. PREPARE DATA FOR GLM
# ============================================================

cat(
  "Preparing data for GLM...\n"
)


# Herbarium records with complete predictor data

datos_model <- datosA_sf[
  !is.na(datosA_sf$dist_road) &
    !is.na(datosA_sf$elevation),
]


# Random points with complete predictor data

random_model <- random_points_sf[
  !is.na(random_points_sf$dist_road) &
    !is.na(random_points_sf$elevation),
]


cat(
  "Records for modelling:",
  nrow(datos_model),
  "\n"
)

cat(
  "Random points for modelling:",
  nrow(random_model),
  "\n"
)


if (
  nrow(datos_model) == 0 ||
  nrow(random_model) == 0
) {
  stop(
    "ERROR: Insufficient data for GLM."
  )
}


# ============================================================
# 17. BINOMIAL GENERALIZED LINEAR MODEL
# ============================================================

cat(
  "Running binomial GLM...\n"
)


data_model <- data.frame(
  tipo = c(
    rep(
      1,
      nrow(datos_model)
    ),

    rep(
      0,
      nrow(random_model)
    )
  ),

  dist = c(
    datos_model$dist_road,
    random_model$dist_road
  ),

  elev = c(
    datos_model$elevation,
    random_model$elevation
  )
)


# ------------------------------------------------------------
# Scale predictors
#
# Scaling allows coefficients to be interpreted relative to
# one standard deviation of each predictor.
# ------------------------------------------------------------

data_model$dist_scaled <- scale(
  data_model$dist
)

data_model$elev_scaled <- scale(
  data_model$elev
)


# ------------------------------------------------------------
# Binomial GLM
# ------------------------------------------------------------

model <- glm(
  tipo ~ dist_scaled + elev_scaled,
  family = "binomial",
  data = data_model,
  na.action = na.omit
)


summary(
  model
)


# ============================================================
# 18. MODEL DIAGNOSTICS
# ============================================================

cat(
  "\n=== MODEL DIAGNOSTICS ===\n"
)


# ------------------------------------------------------------
# Confidence intervals
# ------------------------------------------------------------

conf_int <- confint(
  model
)


cat(
  "\nConfidence intervals:\n"
)

print(
  conf_int
)


# ------------------------------------------------------------
# McFadden's pseudo R-squared
# ------------------------------------------------------------

null_model <- glm(
  tipo ~ 1,
  family = "binomial",
  data = data_model
)


pseudo_r2 <- 1 -
  (
    as.numeric(logLik(model)) /
      as.numeric(logLik(null_model))
  )


cat(
  "\nMcFadden's Pseudo R-squared:",
  pseudo_r2,
  "\n"
)


# ============================================================
# 19. FINAL VALIDATION
# ============================================================

cat(
  "\n=== FINAL DATA SUMMARY ===\n"
)


cat(
  "Original records:",
  nrow(datosA),
  "\n"
)

cat(
  "Records after coordinate cleaning:",
  nrow(datosA_clean),
  "\n"
)

cat(
  "Records within study region:",
  nrow(datosA_sf),
  "\n"
)

cat(
  "Random reference points:",
  nrow(random_points_sf),
  "\n"
)

cat(
  "Records used in GLM:",
  nrow(datos_model),
  "\n"
)

cat(
  "Random points used in GLM:",
  nrow(random_model),
  "\n"
)

# ============================================================
# 20. SESSION INFORMATION
# ============================================================

cat(
  "\n=== SESSION INFORMATION ===\n"
)


sessionInfo()


cat(
  "\nAnalysis complete.\n"
)
