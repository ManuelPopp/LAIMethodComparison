# Packages----------------------------------------------------------------------
require("terra")
require("tidyterra")
require("dplyr")
require("readxl")
require("ggplot2")
library("foreach")
library("doParallel")


# Settings----------------------------------------------------------------------
file_base <- file.path("/dme", "eunomia", "plotdata")

m3m_rasters <- c(
  "gndvi", "green", "gsddsm", "lci", "ndre", "ndvi", "nir", "osavi", "red",
  "rededge", "rgb"
  )

m400_rasters <- c("dom", "dsm", "dtm")


# Functions---------------------------------------------------------------------
## Load spatial data
get_plot_boundary <- function(plot_id) {
  plot_id <- tolower(plot_id)
  
  f_boundary <- file.path(
    file_base, "metadata", "boundaries",
    paste0(plot_id, "_boundary.gpkg")
    )
  
  return(terra::vect(f_boundary))
}


get_raster <- function(plot_id, raster_type) {
  plot_id <- tolower(plot_id)
  raster_type <- tolower(raster_type)
  
  if (raster_type %in% m3m_rasters) {
    sensor_type <- "multispectral"
  } else if (raster_type %in% m400_rasters) {
    sensor_type <- "lidar"
  }
  
  f_raster <- file.path(
    file_base, "measurments", "airborne", sensor_type, raster_type, plot_id,
    paste0(plot_id, "_", raster_type, ".tif")
    )
  
  return(terra::rast(f_raster))
}


## Load LAI data----------------------------------------------------------------
get_lai <- function(plot_id) {
  plot_id <- tolower(plot_id)
  
  f_lai <- file.path(
    file_base, "measurments", "ground", "par", "licor", plot_id,
    paste0(plot_id, ".csv")
  )
  
  return(read.csv(f_lai))
}


## Extract values---------------------------------------------------------------
extract_values <- function(raster, locations, radius = 3) {
  circles <- terra::buffer(locations, width = radius) %>%
    terra::project(raster)
  
  extracted <- terra::extract(
    raster, circles, fun = "mean", cells = FALSE, ID = FALSE, na.rm = TRUE
    )
  
  return(extracted)
}


extract_raster_mean <- function(raster, boundary) {
  bound <- terra::project(boundary, raster)
  
  extracted <- terra::extract(
    raster, bound, fun = "mean", cells = FALSE, ID = FALSE, na.rm = TRUE
  )
  
  return(extracted)
}


## Helper functions-------------------------------------------------------------
cor.coef <- function(x, y) {
  coef <- cor.test(x, y, method = "pearson")$estimate
  return(coef)
}


p.value = function(x, y) {
  p <- cor.test(x, y, method = "pearson")$p.value
  return(p)
}
