# Spatial
library(terra)
library(sf)
library(exactextractr)

# Terrain analysis
library(whitebox)   # for VRM, flow accumulation
library(spatialEco)

# Diversity
library(vegan)
library(betapart)

# Modelling
library(mgcv)
library(gratia)

# Data handling
library(tidyverse)

# plots with lat-long csv
plots_df <- read.csv("Data/plot_coordinates.csv")
head(plots_df)

plots <- st_as_sf(
  plots_df,
  coords = c("longitude", "latitude"),
  crs = 4326
)

dem <- rast("DEM/ALOS_PALSAR_12p5m_DEM.tif")

# Transform plot points to DEM CRS
plots <- st_transform(plots, crs(dem))

plots_buf <- st_buffer(plots, dist = 56.42)
# Basic terrain derivatives
slope  <- terrain(dem, v = "slope", unit = "degrees")
aspect <- terrain(dem, v = "aspect", unit = "degrees")
tri    <- terrain(dem, v = "TRI")
rough  <- terrain(dem, v = "roughness")
#Topographic Position Index (TPI)
# dissolve all square plots
plots_union <- st_union(plots_buf)

# convert to terra vector
plots_vect <- vect(plots_union)

# crop DEM
dem_crop <- crop(dem, plots_vect)
dem_crop <- mask(dem_crop, plots_vect)

tpi_25 <- tpi(dem_crop, scale = 25)

# tpi <- tpi(dem, scale = 25)   # ~300 m window (appropriate for Ladakh)


# Hydrological variables
dem_filled <- fill(dem)

flow_acc <- flowAccumulation(dem_filled)
twi <- terrain(dem_filled, v = "TWI")

# Vector Ruggedness Measure (VRM)
writeRaster(dem, "dem.tif", overwrite = TRUE)

wbt_vector_ruggedness_measure(
  dem = "dem.tif",
  output = "vrm.tif",
  window = 3
)

vrm <- rast("vrm.tif")

# STACK ALL MICROTOPOGRAPHIC VARIABLES
topo_stack <- c(
  dem,
  slope,
  aspect,
  tri,
  rough,
  tpi,
  twi,
  flow_acc,
  vrm
)

names(topo_stack) <- c(
  "elevation",
  "slope",
  "aspect",
  "tri",
  "roughness",
  "tpi",
  "twi",
  "flow_acc",
  "vrm"
)

# EXTRACT MEAN + SD WITHIN EACH 1-HA PLOT

# This is where microhabitat heterogeneity enters

micro_mean <- exact_extract(topo_stack, plots_buf, "mean")
micro_sd   <- exact_extract(topo_stack, plots_buf, "stdev")

# FORMAT MICROTOPOGRAPHY TABLE
micro <- data.frame(
  Plot_ID = plots$Plot_ID,
  micro_mean,
  micro_sd
)

# Rename columns clearly
names(micro) <- c(
  "Plot_ID",
  paste0(names(topo_stack), "_mean"),
  paste0(names(topo_stack), "_sd")
)


# (Optional CV for selected variables)

micro <- micro %>%
  mutate(
    slope_cv = slope_sd / slope_mean,
    twi_cv   = twi_sd / twi_mean,
    tpi_cv   = tpi_sd / abs(tpi_mean)
  )

#   SPECIES ABUNDANCE MATRIX
comm <- read.csv("species_abundance.csv", row.names = 1)

# HILL DIVERSITY (ALPHA)
hill_div <- data.frame(
  Plot_ID = rownames(comm),
  Hill_q0 = specnumber(comm),
  Hill_q1 = exp(diversity(comm, index = "shannon")),
  Hill_q2 = 1 / diversity(comm, index = "simpson"))

# BETA DIVERSITY (SORRENSEN)
beta <- beta.pair(comm, index.family = "sorensen")

beta_df <- data.frame(
  Plot_ID = rownames(comm),
  beta_turnover = rowMeans(as.matrix(beta$beta.sim), na.rm = TRUE),
  beta_nestedness = rowMeans(as.matrix(beta$beta.sne), na.rm = TRUE),
  beta_total = rowMeans(as.matrix(beta$beta.sor), na.rm = TRUE)
)


# MERGE ALL DATA
div_data <- micro %>%
  left_join(hill_div, by = "Plot_ID") %>%
  left_join(beta_df, by = "Plot_ID")

# GAMs WITH ECOLOGICALLY SUITABLE FAMILIES
Hill q0 (richness → counts)
gam_q0 <- gam(
  Hill_q0 ~
    s(elevation_mean, k = 5) +
    s(tpi_sd, k = 5) +
    s(twi_sd, k = 5) +
    s(roughness_sd, k = 5),
  family = nb(),
  data = div_data,
  method = "REML"
)

# Hill q1 (effective Shannon → continuous)
gam_q1 <- gam(
  Hill_q1 ~
    s(tpi_sd, k = 5) +
    s(twi_sd, k = 5) +
    s(slope_sd, k = 5) +
    s(vrm_sd, k = 5),
  family = Gamma(link = "log"),
  data = div_data,
  method = "REML"
)

# Beta turnover (bounded 0–1)
div_data$beta_turnover <-
  pmin(pmax(div_data$beta_turnover, 0.001), 0.999)

gam_beta <- gam(
  beta_turnover ~
    s(elevation_mean, k = 5) +
    s(tpi_sd, k = 5) +
    s(twi_sd, k = 5),
  family = betar(link = "logit"),
  data = div_data,
  method = "REML"
)

#  VISUALISATION
draw(gam_q1, residuals = TRUE, rug = TRUE)
draw(gam_beta, residuals = TRUE, rug = TRUE)
