# Load required libraries
library("sf")
library("terra")
library("tidyterra")
library("lidR")

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
recompute <- "--recompute" %in% args
plot_args <- args[args != "--recompute"]

if (length(plot_args) != 1) {
  print(
    "Plot ID is a required argument which can be either an index or a plot ID."
    )
  stop("Usage: Rscript lidar_metrics.R <plot_id>")
}

# Input directories
dir_base <- "/dme/eunomia/plotdata"
dir_meas <- file.path(dir_base, "measurments")
dir_meta <- file.path(dir_base, "metadata")
dir_air <- file.path(dir_meas, "airborne")
dir_lidar <- file.path(dir_air, "lidar")
dir_las <- file.path(dir_lidar, "pointclouds")
dir_dtm <- file.path(dir_lidar, "dem")
dir_dsm <- file.path(dir_lidar, "dsm")
dir_bnds <- file.path(dir_meta, "boundaries_adjusted")

if (grepl("^[+-]?[0-9]+$", plot_args[1])) {
  plot_list <- list.dirs(dir_lidar, full.names = FALSE)
  plot_id <- plot_list[as.numeric(plot_args[1])]
} else {
  plot_id <- plot_args[1]
}

src_las <- file.path(dir_las, plot_id, paste0(plot_id, "_pointcloud.las"))
src_bnd <- file.path(dir_bnds, paste0(plot_id, "_adjusted.gpkg"))
src_dtm <- file.path(dir_dtm, plot_id, paste0(plot_id, "_dem.tif"))
src_dsm <- file.path(dir_dsm, plot_id, paste0(plot_id, "_dsm.tif"))

# Output directories
dir_dst <- "/storage/poppman/lai_metrics"
dir_rst <- file.path(dir_dst, "metric_rasters")
dir_csv <- file.path(dir_dst, "metric_tables")

dst_csv <- file.path(dir_csv, paste0(plot_id, "_lidar_metrics.csv"))

# Stop if metrics have already been computed and recompute is not requested
if (file.exists(dst_csv) && !recompute) {
  cat("Metrics already computed for plot:", plot_id, "\n")
  quit(save = "no", status = 0)
}

# Height threshold for tree points (in meters)
th_z <- 1.5

# Load input data
if (!file.exists(src_las)) {
  stop(paste("The LAS file does not exist:", src_las))
}
if (!file.exists(src_bnd)) {
  stop(paste("The boundary file does not exist:", src_bnd))
}
if (!file.exists(src_dtm)) {
  stop(paste("The DTM file does not exist:", src_dtm))
}
if (!file.exists(src_dsm)) {
  stop(paste("The DSM file does not exist:", src_dsm))
}

cat("Processing plot:", plot_id, "\n")

las <- lidR::readALS(src_las)
boundary <- sf::st_read(src_bnd, quiet = TRUE)
dtm <- terra::rast(src_dtm)
dsm <- terra::rast(src_dsm)

if (is.empty(las)) {
  stop("The LAS file contains no points.")
}

# Get CRS from LAS and transform boundary if necessary
crs_las <- lidR::st_crs(las)
crs_bnd <- sf::st_crs(boundary)

if (crs_las != crs_bnd) {
  boundary <- sf::st_transform(boundary, crs_las)
}

# Clip and normalise the point cloud
las_clipped <- lidR::clip_roi(las, boundary)

if (is.empty(las_clipped)) {
  stop("The clipped point cloud contains no points.")
}

las_normalised <- lidR::normalize_height(las_clipped, dtm)

if (is.empty(las_normalised)) {
  stop("The normalised point cloud contains no points.")
}

# Compute the above and below ratio index (ABRI) sumnall2016
print("Computing ABRI...")
dst_abri <- file.path(dir_rst, "abri", paste0(plot_id, "_abri.tif"))
if (!file.exists(dst_abri) || recompute) {
  dir.create(dir_abri, recursive = TRUE, showWarnings = FALSE)

  lidR::grid_metrics(
    las_normalised,
    ~ sum(Z >= th_z) / sum(Z < th_z),
    res = 1
  ) %>%
    terra::writeRaster(
      dst_abri,
      overwrite = TRUE
    )
}

abri <- terra::rast(dst_abri)
abri_mean <- terra::global(abri, fun = "mean", na.rm = TRUE)
abri_median <- terra::global(abri, fun = median, na.rm = TRUE)
abri_sd <- terra::global(abri, fun = "sd", na.rm = TRUE)
rm(abri)

# Compute the laser penetration index 1
print("Computing LP1...")
dst_lp1 <- file.path(dir_rst, "lp1", paste0(plot_id, "_lp1.tif"))
if (!file.exists(dst_lp1) || recompute) {
  dir.create(dir_lp1, recursive = TRUE, showWarnings = FALSE)

  lidR::grid_metrics(
    las_normalised,
    ~ sum(Z < th_z) / length(Z),
    res = 1
  ) %>%
    terra::writeRaster(
      dst_lp1,
      overwrite = TRUE
    )
}

lp1 <- terra::rast(dst_lp1)
lp1_mean <- terra::global(lp1, fun = "mean", na.rm = TRUE)
lp1_median <- terra::global(lp1, fun = median, na.rm = TRUE)
lp1_sd <- terra::global(lp1, fun = "sd", na.rm = TRUE)
rm(lp1)

# Compute the laser penetration index 2
print("Computing LP2...")
lp2 <- lidR::grid_metrics(
  las_normalised,
  ~ (
    sum(Z < th_z & ReturnNumber == 1) + sum(Z < th_z & ReturnNumber == NumberOfReturns)
    ) / (
        0.5 * (sum(ReturnNumber == 1) + sum(ReturnNumber == NumberOfReturns))
        ),
  res = 1
)

dst_lp2 <- file.path(dir_rst, "lp2", paste0(plot_id, "_lp2.tif"))
if (!file.exists(dst_lp2) || recompute) {
  dir.create(dir_lp2, recursive = TRUE, showWarnings = FALSE)

  lidR::grid_metrics(
    las_normalised,
    ~ (
      sum(Z < th_z & ReturnNumber == 1) + sum(Z < th_z & ReturnNumber == NumberOfReturns)
      ) / (
          0.5 * (sum(ReturnNumber == 1) + sum(ReturnNumber == NumberOfReturns))
          ),
    res = 1
  ) %>%
    terra::writeRaster(
      dst_lp2,
      overwrite = TRUE
    )
}

lp2 <- terra::rast(dst_lp2)
lp2_mean <- terra::global(lp2, fun = "mean", na.rm = TRUE)
lp2_median <- terra::global(lp2, fun = median, na.rm = TRUE)
lp2_sd <- terra::global(lp2, fun = "sd", na.rm = TRUE)
rm(lp2)

# Compute hanopy height model
print("Computing CHM...")
chm <- (dsm - terra::resample(dtm, dsm)) %>%
  terra::crop(terra::project(terra::vect(boundary), dsm))

dst_chm <- file.path(dir_rst, "chm", paste0(plot_id, "_chm.tif"))
if (!file.exists(dst_chm) || recompute) {
  dir.create(dst_chm, recursive = TRUE, showWarnings = FALSE)

  (dsm - terra::resample(dtm, dsm)) %>%
    terra::crop(terra::project(terra::vect(boundary), dsm)) %>%
    terra::writeRaster(
      dst_chm,
      overwrite = TRUE
    )
}

chm <- terra::rast(dst_chm)
chm_mean <- terra::global(chm, fun = "mean", na.rm = TRUE)
chm_median <- terra::global(chm, fun = median, na.rm = TRUE)
chm_sd <- terra::global(chm, fun = "sd", na.rm = TRUE)
rm(chm)

# Save metrics to CSV
df_metrics <- data.frame(
  plot_id = plot_id,
  abri_mean = abri_mean,
  abri_median = abri_median,
  abri_sd = abri_sd,
  lp1_mean = lp1_mean,
  lp1_median = lp1_median,
  lp1_sd = lp1_sd,
  lp2_mean = lp2_mean,
  lp2_median = lp2_median,
  lp2_sd = lp2_sd,
  chm_mean = chm_mean,
  chm_median = chm_median,
  chm_sd = chm_sd
)

write.csv(df_metrics, dst_csv, row.names = FALSE)

print("Finished.")