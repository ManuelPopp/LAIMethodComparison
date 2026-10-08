library("sf")
library("terra")
library("tidyterra")
library("lidR")

src_las <- "/dme/eunomia/plotdata/measurments/airborne/lidar/pointclouds/bl01/bl01_pointcloud.las"
src_bnd <- "/dme/eunomia/plotdata/metadata/boundaries_adjusted/bl01_adjusted.gpkg"
src_dtm <- "/dme/eunomia/plotdata/measurments/airborne/lidar/dem/bl01/bl01_dem.tif"
src_dsm <- "/dme/eunomia/plotdata/measurments/airborne/lidar/dsm/bl01/bl01_dsm.tif"

# Height threshold for tree points (in meters)
th_z <- 1.5

# Load input data
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

if (is.empty(las_normalised)) {
  stop("The clipped point cloud contains no points.")
}

las_normalised <- lidR::normalize_height(las_clipped, dtm)

if (is.empty(las_normalised)) {
  stop("The normalised point cloud contains no points.")
}

# Compute the above and below ratio index (ABRI) sumnall2016
abri <- lidR::grid_metrics(
  las_normalised,
  ~ sum(Z >= th_z) / sum(Z < th_z),
  res = 1
)

# Compute the laser penetration index 1
lp1 <- lidR::grid_metrics(
  las_normalised,
  ~ sum(Z < th_z) / length(Z),
  res = 1
)

# Compute the laser penetration index 2
lp2 <- lidR::grid_metrics(
  las_normalised,
  ~ (
    sum(Z < th_z & ReturnNumber == 1) + sum(Z < th_z & ReturnNumber == NumberOfReturns)
    ) / (
        0.5 * (sum(ReturnNumber == 1) + sum(ReturnNumber == NumberOfReturns))
        ),
  res = 1
)

# Compute hanopy height model
chm <- (dsm - terra::resample(dtm, dsm)) %>%
  terra::crop(terra::project(terra::vect(boundary), dsm))