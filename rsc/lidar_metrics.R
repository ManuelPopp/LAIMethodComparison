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
  plot_list <- sub(
    "_adjusted.gpkg", "",
    list.files(dir_bnds, full.names = FALSE, recursive = FALSE)
  )
  plot_list <- plot_list[plot_list != "."]
  plot_id <- plot_list[as.numeric(plot_args[1])]
} else {
  plot_id <- plot_args[1]
}
print(paste("Processing plot:", plot_id))

src_las <- file.path(dir_las, plot_id, paste0(plot_id, "_pointcloud.las"))
src_bnd <- file.path(dir_bnds, paste0(plot_id, "_adjusted.gpkg"))
src_dtm <- file.path(dir_dtm, plot_id, paste0(plot_id, "_dem.tif"))
src_dsm <- file.path(dir_dsm, plot_id, paste0(plot_id, "_dsm.tif"))

src_dsm_m3m <- file.path(
  sub("lidar", "multispectral", dir_dsm),
  plot_id, paste0(plot_id, "_m3m_dsm.tif")
  )
src_dsm_m4t <- file.path(
  sub("lidar", "thermal", dir_dsm),
  plot_id, paste0(plot_id, "_m4t_dsm.tif")
  )

# Output directories
dir_dst <- "/storage/poppman/lai_metrics"
dir_rst <- file.path(dir_dst, "metric_rasters")
dir_csv <- file.path(dir_dst, "metric_tables")
dir_data <- file.path(dir_dst, "intermediate_data", plot_id)

dst_csv <- file.path(dir_csv, paste0(plot_id, "_lidar_metrics.csv"))

dir.create(dir_csv, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_rst, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_data, recursive = TRUE, showWarnings = FALSE)

# Stop if metrics have already been computed and recompute is not requested
if (file.exists(dst_csv) && !recompute) {
  cat("Metrics already computed for plot:", plot_id, "\n")
  quit(save = "no", status = 0)
}

# Height threshold for tree points (in meters) and tree cover
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

dsm_l2 <- terra::rast(src_dsm)
dtm <- terra::rast(src_dtm) %>%
  terra::project(dsm_l2, method = "bilinear")

boundary <- sf::st_read(src_bnd, quiet = TRUE) %>%
  sf::st_transform(crs = terra::crs(dsm_l2))

if (file.exists(src_dsm_m3m)) {
  dsm_m3m <- terra::rast(src_dsm_m3m) %>%
    terra::project(dsm_l2, method = "bilinear")
} else {
  warning(paste("The DSM M3M file does not exist:", src_dsm_m3m))
  dsm_m3m <- NULL
}
if (file.exists(src_dsm_m4t)) {
  dsm_m4t <- terra::rast(src_dsm_m4t) %>%
  terra::project(dsm_l2, method = "bilinear")
} else {
  warning(paste("The DSM M4T file does not exist:", src_dsm_m4t))
  dsm_m4t <- NULL
}

#>---------------------------------------------------------------------<
#> Functions
#>---------------------------------------------------------------------<
# Replace invalid raster values (Inf, -Inf, NaN) with NA
replace_nonfinite <- function(r) {
  if (!inherits(r, "SpatRaster")) {
    r <- terra::rast(r)
  }
  terra::app(
    r,
    fun = function(x) {
      x[!is.finite(x)] <- NA_real_
      x
    }
  )
}

#>---------------------------------------------------------------------<
#> Prepare data
#>---------------------------------------------------------------------<
# LiDAR point cloud
## Function to load point cloud when required
get_las <- function(las_file = src_las) {
  las <- lidR::readALS(las_file)
  las_filtered <- lidR::filter_duplicates(las)
  rm(las)
  if (is.empty(las_filtered)) {
    stop("The LAS file contains no points.")
  }

  # Get CRS from LAS and transform boundary if necessary
  crs_las <- lidR::st_crs(las_filtered)
  crs_bnd <- sf::st_crs(boundary)

  # Clip and normalise the point cloud
  if (!sf::st_crs(crs_las) == sf::st_crs(crs_bnd)) {
    warning("CRS of LAS and boundary do not match.")
    boundary <- sf::st_transform(boundary, crs = crs_las)
  }
  las_clipped <- lidR::clip_roi(
    las_filtered, sf::st_buffer(boundary, dist = 1)
    )

  if (is.empty(las_clipped)) {
    stop("The clipped point cloud contains no points.")
  }
  las_clipped <<- las_clipped
  las_normalised <<- lidR::normalize_height(las_clipped, dtm)

  if (is.empty(las_normalised)) {
    stop("The normalised point cloud contains no points.")
  }
}

# Raster data
## Prepare DSM by filling missing values and smoothing
if (is.null(dsm_m3m) && is.null(dsm_m4t)) {
  stop("Both DSM M3M and DSM M4T files are missing. Cannot proceed.")
}

dir.create(dir_data, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_rst, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(file.path(dir_data, "dsm_p2r.tif")) || recompute) {
  get_las()
  lidR::rasterize_canopy(
    las_clipped,
    res = dsm_l2,
    algorithm = lidR::p2r()
  ) %>%
    terra::writeRaster(
      filename = file.path(dir_data, "dsm_p2r.tif"),
      overwrite = TRUE
    )
  cat(
    "\nPoint cloud rasterised. Wrote to:",
    file.path(dir_data, "dsm_p2r.tif")
    )
}

if (!file.exists(file.path(dir_data, "dsm_filled.tif")) || recompute) {
  cat("\nFilling missing values in DSM...")
  if (!is.null(dsm_m3m) && !is.null(dsm_m4t)) {
    gap_fill <- median(dsm_m3m, dsm_m4t, dsm_l2, na.rm = TRUE)
  } else if (!is.null(dsm_m4t)) {
    gap_fill <- median(dsm_m4t, dsm_l2, na.rm = TRUE)
  } else if (!is.null(dsm_m3m)) {
    gap_fill <- median(dsm_m3m, dsm_l2, na.rm = TRUE)
  } else {
    stop("Both DSM M3M and DSM M4T files are missing. Cannot fill gaps.")
  }

  terra::rast(file.path(dir_data, "dsm_p2r.tif")) %>%
    terra::cover(
      gap_fill,
      filename = file.path(dir_data, "dsm_filled.tif"),
      overwrite = TRUE
    )
  cat(
    "\nDSM filling completed. Wrote to:",
    file.path(dir_data, "dsm_filled.tif")
    )
}

if (!file.exists(file.path(dir_data, "dsm_smooth.tif")) || recompute) {
  cat("\nSmoothing DSM...")
  terra::focal(
    terra::rast(file.path(dir_data, "dsm_filled.tif")),
    w = 3,
    fun = "mean",
    na.rm = TRUE,
    filename = file.path(dir_data, "dsm_smooth.tif"),
    overwrite = TRUE
  )
  cat(
    "\nDSM smoothing completed. Wrote to:",
    file.path(dir_data, "dsm_smooth.tif")
    )
}

# Load final DSM
dsm <- terra::rast(file.path(dir_data, "dsm_smooth.tif"))

#>---------------------------------------------------------------------<
#> Compute metrics
#>---------------------------------------------------------------------<
# Compute the above and below ratio index (ABRI) sumnall2016
cat("\nComputing ABRI...")
dst_abri <- file.path(dir_rst, "abri", paste0(plot_id, "_abri.tif"))
if (!file.exists(dst_abri) || recompute) {
  dir.create(
    file.path(dir_rst, "abri"), recursive = TRUE, showWarnings = FALSE
    )
  get_las()
  lidR::grid_metrics(
    las_normalised,
    ~ sum(Z >= th_z) / sum(Z < th_z),
    res = 1
  ) %>%
  replace_nonfinite() %>%
    terra::writeRaster(
      filename = dst_abri,
      overwrite = TRUE
    )
}

abri <- terra::rast(dst_abri) %>%
  terra::crop(terra::project(terra::vect(boundary), .))

abri_mean <- terra::global(abri, fun = "mean", na.rm = TRUE)[1, 1]
abri_median <- terra::global(abri, fun = median, na.rm = TRUE)[1, 1]
abri_sd <- terra::global(abri, fun = "sd", na.rm = TRUE)[1, 1]
rm(abri)

# Compute the laser penetration index 1
cat("\nComputing LP1...")
dst_lp1 <- file.path(dir_rst, "lp1", paste0(plot_id, "_lp1.tif"))
if (!file.exists(dst_lp1) || recompute) {
  dir.create(file.path(dir_rst, "lp1"), recursive = TRUE, showWarnings = FALSE)
  get_las()
  lidR::grid_metrics(
    las_normalised,
    ~ sum(Z < th_z) / length(Z),
    res = 1
  ) %>%
    replace_nonfinite() %>%
    terra::writeRaster(
      filename = dst_lp1,
      overwrite = TRUE
    )
}

lp1 <- terra::rast(dst_lp1) %>%
  terra::crop(terra::project(terra::vect(boundary), .))

lp1_mean <- terra::global(lp1, fun = "mean", na.rm = TRUE)[1, 1]
lp1_median <- terra::global(lp1, fun = median, na.rm = TRUE)[1, 1]
lp1_sd <- terra::global(lp1, fun = "sd", na.rm = TRUE)[1, 1]
rm(lp1)

# Compute the laser penetration index 2
cat("\nComputing LP2...")
dst_lp2 <- file.path(dir_rst, "lp2", paste0(plot_id, "_lp2.tif"))
if (!file.exists(dst_lp2) || recompute) {
  dir.create(file.path(dir_rst, "lp2"), recursive = TRUE, showWarnings = FALSE)
  get_las()
  lidR::grid_metrics(
    las_normalised,
    ~ (
      sum(Z < th_z & ReturnNumber == 1) +
      sum(Z < th_z & ReturnNumber == NumberOfReturns)
      ) / (
          0.5 * (sum(ReturnNumber == 1) + sum(ReturnNumber == NumberOfReturns))
          ),
    res = 1
  ) %>%
  replace_nonfinite() %>%
    terra::writeRaster(
      filename = dst_lp2,
      overwrite = TRUE
    )
}

lp2 <- terra::rast(dst_lp2) %>%
  terra::crop(terra::project(terra::vect(boundary), .))

lp2_mean <- terra::global(lp2, fun = "mean", na.rm = TRUE)[1, 1]
lp2_median <- terra::global(lp2, fun = median, na.rm = TRUE)[1, 1]
lp2_sd <- terra::global(lp2, fun = "sd", na.rm = TRUE)[1, 1]
rm(lp2)

# Compute the intensity ratio index (IRI)
cat("\nComputing IRI...")
dst_iri <- file.path(dir_rst, "iri", paste0(plot_id, "_iri.tif"))
if (!file.exists(dst_iri) || recompute) {
  dir.create(file.path(dir_rst, "iri"), recursive = TRUE, showWarnings = FALSE)
  get_las()
  lidR::grid_metrics(
    las_normalised,
    ~ mean(Intensity[Z >= th_z], na.rm = TRUE) /
      mean(Intensity[Z < th_z], na.rm = TRUE),
    res = 1
  ) %>%
  replace_nonfinite() %>%
    terra::writeRaster(
      filename = dst_iri,
      overwrite = TRUE
    )
}

# Compute hanopy height model
cat("\nComputing CHM...")
dst_chm <- file.path(dir_rst, "chm", paste0(plot_id, "_chm.tif"))
if (!file.exists(dst_chm) || recompute) {
  dir.create(file.path(dir_rst, "chm"), recursive = TRUE, showWarnings = FALSE)
  (dsm - dtm) %>%
    terra::crop(terra::project(terra::vect(boundary), dsm)) %>%
    terra::writeRaster(
      filename = dst_chm,
      overwrite = TRUE
    )
}

chm <- terra::rast(dst_chm) %>%
  terra::crop(terra::project(terra::vect(boundary), .))

chm_mean <- terra::global(chm, fun = "mean", na.rm = TRUE)[1, 1]
chm_median <- terra::global(chm, fun = median, na.rm = TRUE)[1, 1]
chm_sd <- terra::global(chm, fun = "sd", na.rm = TRUE)[1, 1]
chm_iqr <- terra::global(chm, fun = IQR, na.rm = TRUE)[1, 1]
chm_p90 <- terra::global(
  chm,
  fun = function(x, ...) {quantile(x, probs = 0.90, na.rm = TRUE)}
  )[1, 1]

# Tree cover
dst_treecover <- file.path(
  dir_rst, "treecover", paste0(plot_id, "_treecover.tif")
  )

if (!file.exists(dst_treecover) || recompute) {
  dir.create(
    file.path(dir_rst, "treecover"), recursive = TRUE, showWarnings = FALSE
    )
  (chm >= th_z) %>%
    terra::writeRaster(
      filename = dst_treecover,
      overwrite = TRUE
    )
}

tree_cover <- terra::rast(dst_treecover) %>%
  terra::crop(terra::project(terra::vect(boundary), .))

tree_cover_fraction <- terra::global(
  tree_cover, fun = "mean", na.rm = TRUE
  )[1, 1]


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
  chm_sd = chm_sd,
  chm_p90 = chm_p90,
  chm_iqr = chm_iqr,
  tree_cover_fraction = tree_cover_fraction
)

write.csv(df_metrics, dst_csv, row.names = FALSE)

cat("\nFinished.")