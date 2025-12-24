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
plots_df <- read.csv("plot_coordinates.csv")
head(plots_df)

plots <- st_as_sf(
  plots_df,
  coords = c("longitude", "latitude"),
  crs = 4326
)
