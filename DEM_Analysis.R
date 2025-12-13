##############################################################
### 1. LOAD REQUIRED PACKAGES
##############################################################
packages <- c("terra","sf","dplyr","vegan","betapart",
              "mgcv","gratia","car","spatialEco","insol","ggplot2")
installed <- rownames(installed.packages())
for (p in packages) if (!(p %in% installed)) install.packages(p)
lapply(packages, library, character.only = TRUE)
install.packages("remotes")
remotes::install_version(
  "insol",
  version = "1.2.1",
  repos = "https://cran.r-project.org"
)

##############################################################
### 2. IMPORT SPECIES × SITE MATRIX AND COORDINATES
##############################################################
species_data <- read.csv("Data/species_site_matrix.csv", row.names = 1)
coords <- read.csv("Data/plot_coordinates.csv")  # Must have Plot_ID, Latitude, Longitude

##############################################################
### 3. CREATE SPATIAL OBJECTS & REPROJECT TO UTM
##############################################################
library(sf)

plots_sf <- st_as_sf(coords, coords = c("Longitude", "Latitude"), crs = 4326)
utm_zone <- floor((mean(coords$Longitude) + 180)/6) + 1
utm_crs <- paste0("+proj=utm +zone=", utm_zone, " +datum=WGS84 +units=m +no_defs")
plots_sf <- st_transform(plots_sf, crs = utm_crs)

##############################################################
### 4. LOAD DEM AND DERIVE MICROTOPOGRAPHIC VARIABLES
##############################################################
library(terra)

dem <- rast("DEM/ALOS_PALSAR_12p5m_DEM.tif")
dem <- project(dem, utm_crs)

# Primary terrain metrics
slope <- terra::terrain(dem, v = "slope", unit = "degrees")
aspect  <- terra::terrain(dem, v = "aspect", unit = "degrees")

# Secondary terrain metrics
#terra::terrain
#profile_curv <- terra::terrain(dem, v = "profilecurvature")
#plan_curv    <- terra::terrain(dem, v = "plancurvature")
#tpi          <- terra::terrain(dem, v = "tpi")
#tri          <- terra::terrain(dem, v = "tri")

roughness    <- terra::terrain(dem, v = "roughness")
library(spatialEco)

tpi <- spatialEco::tpi(dem)
tri <- spatialEco::tri(dem)

# Hydrological metric
flowdir <- terra::terrain(dem, v = "flowdir")
flowacc <- terra::flowAccumulation(flowdir)
flowacc_log <- log1p(flowacc)
slope_rad <- slope * pi / 180.0
twi <- log((flowacc + 1) / tan(slope_rad))
names(twi) <- "twi"

# Ruggedness
vrm <- spatialEco::vrm(dem, s = 3)
vrm <- rast(vrm); names(vrm) <- "vrm"

##############################################################
### 5. MONTHLY SOLAR RADIATION (JUNE–SEPT)
##############################################################
#slope_rad <- radians(values(slope))
aspect_rad <- values(aspect)*pi/180.0
lat <- mean(coords$Latitude)
lat_rad <- lat* pi/180.0
library(insol)
get_solar_month <- function(DOY){
  sol <- insolation(slope = slope_rad, aspect = aspect_rad, lat = lat_rad,
                    J = DOY, local = TRUE, hourly = FALSE)
  r <- rast(dem); values(r) <- sol; return(r)
}
days <- c(June=172, July=203, August=234, September=265)
solar_june <- get_solar_month(days["June"])
solar_july <- get_solar_month(days["July"])
solar_aug  <- get_solar_month(days["August"])
solar_sep  <- get_solar_month(days["September"])

solar_stack <- c(solar_june, solar_july, solar_aug, solar_sep)
names(solar_stack) <- c("solar_June","solar_July","solar_August","solar_September")

##############################################################
### 6. STACK VARIABLES
##############################################################
#micro_stack <- c(dem, slope, aspect, tpi, tri, roughness, twi, vrm, solar_stack)
#names(micro_stack)[1:8] <- c("elevation","slope","aspect","tpi","tri",
#                             "roughness","twi","vrm")

micro_stack <- c(dem, slope, aspect, tpi, tri, roughness, twi, vrm)

names(micro_stack)[1:8] <- c("elevation","slope","aspect","tpi","tri",
                             "roughness","twi","vrm")
##############################################################
### 7. 1-HA BUFFERS AND EXTRACTION (MEAN, SD, CV)
##############################################################
plots_buffer <- st_buffer(plots_sf, dist = 56.4)

# Function to compute mean, SD, CV
extract_mean_sd_cv <- function(rast_layer, buffers, var_name) 
  {
  df_mean <- terra::extract(rast_layer, buffers, fun = mean, na.rm = TRUE)
  df_sd   <- terra::extract(rast_layer, buffers, fun = sd, na.rm = TRUE)
  df <- data.frame(
    Plot_ID = coords$Plot_ID,
    mean = df_mean[,2],
    sd = df_sd[,2]
  )
  df$cv <- df$sd / df$mean
  names(df)[2:4] <- paste0(var_name, c("_mean", "_sd", "_cv"))
  return(df)
}

# Variables benefiting most from SD/CV
var_list <- list(
  extract_mean_sd_cv(slope, plots_buffer, "slope"),
  extract_mean_sd_cv(twi, plots_buffer, "twi"),
  extract_mean_sd_cv(tri, plots_buffer, "tri"),
  extract_mean_sd_cv(tpi, plots_buffer, "tpi"),
  extract_mean_sd_cv(roughness, plots_buffer, "roughness")
)

topo_het <- Reduce(function(x, y) merge(x, y, by = "Plot_ID"), var_list)

##############################################################
### 8. ALPHA & BETA DIVERSITY
##############################################################
library(vegan)
alpha_div <- data.frame(
  Plot_ID = rownames(species_data),
  Richness = rowSums(species_data > 0),
  Shannon  = diversity(species_data, index = "shannon"),
  Simpson  = diversity(species_data, index = "simpson")
)
library(betapart)
beta_parts <- beta.pair(species_data, index.family = "sorensen")
beta_summary <- data.frame(
  Total_Beta_Sorensen = mean(beta_parts$beta.sor),
  Turnover_Beta_Sim = mean(beta_parts$beta.sim),
  Nestedness_Beta_Sne = mean(beta_parts$beta.sne)
)

div_micro <- merge(alpha_div, topo_het, by = "Plot_ID")

##############################################################
### 9. CHECK MULTICOLLINEARITY
##############################################################
vif_mod <- lm(Shannon ~ slope_mean + slope_sd + twi_mean + twi_sd +
              tri_mean + tri_sd + tpi_mean + tpi_sd + roughness_mean + roughness_sd,
              data = div_micro)
vif_vals <- car::vif(vif_mod)
print(vif_vals)

# Optionally remove variables with VIF > 5
div_micro <- div_micro[, !names(div_micro) %in% names(vif_vals[vif_vals > 5])]

##############################################################
### 10. FIT GAM MODEL (MEAN + SD VARIABLES)
##############################################################
library(mgcv)
gam_shannon <- gam(Shannon ~ 
                     s(slope_mean, k=5) + s(slope_sd, k=5) +
                     s(twi_mean, k=5)  + s(twi_sd, k=5) +
                     s(tri_mean, k=5)  + s(tri_sd, k=5) +
                     s(tpi_mean, k=5)  + s(tpi_sd, k=5) +
                     s(roughness_mean, k=5) + s(roughness_sd, k=5),
                   data = div_micro, method = "REML")

summary(gam_shannon)
gam.check(gam_shannon)
cat("\nR² =", summary(gam_shannon)$r.sq,
    "\nDeviance explained =", round(summary(gam_shannon)$dev.expl * 100, 2), "%\n")

##############################################################
### 11. VARIABLE IMPORTANCE (EDF & F-VALUES)
##############################################################
imp <- summary(gam_shannon)$s.table
importance <- data.frame(
  Variable = rownames(imp),
  EDF = imp[,"edf"],
  F_value = imp[,"F"],
  p_value = imp[,"p-value"]
)
importance <- importance[order(-importance$F_value), ]
write.csv(importance, "GAM_variable_importance.csv", row.names = FALSE)

ggplot(importance, aes(x=reorder(Variable, F_value), y=F_value)) +
  geom_col() + coord_flip() +
  theme_minimal() +
  labs(title="Relative Influence of Predictors (GAM F-values)",
       x="Predictor", y="F-value")

##############################################################
### 12. VISUALIZATION
##############################################################
par(mfrow=c(3,4))
plot(gam_shannon, shade=TRUE, seWithMean=TRUE, scale=0)

draw(gam_shannon, residuals = TRUE)
vis.gam(gam_shannon, view = c("twi_sd", "slope_sd"), plot.type = "contour",
        color = "topo", main = "Interaction: TWI_SD × Slope_SD on Shannon Diversity")

##############################################################
### 13. SAVE OUTPUTS
##############################################################
write.csv(alpha_div, "alpha_diversity_results.csv", row.names=FALSE)
write.csv(beta_summary, "beta_diversity_summary.csv", row.names=FALSE)
write.csv(div_micro, "diversity_microtopography_SD_CV.csv", row.names=FALSE)
saveRDS(gam_shannon, "GAM_Shannon_SD_CV_Model.rds")

##############################################################
### END OF SCRIPT
##############################################################
