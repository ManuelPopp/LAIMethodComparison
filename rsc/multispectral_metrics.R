# Load required libraries
library("sf")
library("terra")
library("tidyterra")

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

# Updatevalue to set non-tree pixels to NA or zero
udv <- 0

# Input directories
dir_base <- "/dme/eunomia/plotdata"
dir_meas <- file.path(dir_base, "measurments")
dir_meta <- file.path(dir_base, "metadata")
dir_air <- file.path(dir_meas, "airborne")
dir_ms <- file.path(dir_air, "multispectral")
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

src_bnd <- file.path(dir_bnds, paste0(plot_id, "_adjusted.gpkg"))

# Intermediate files (output from LiDAR processing)
dir_dst <- "/storage/poppman/lai_metrics"
dir_rst <- file.path(dir_dst, "metric_rasters")
src_treecover <- file.path(
  dir_rst, "treecover", paste0(plot_id, "_treecover.tif")
  )

# Output directories
dir_csv <- file.path(dir_dst, "metric_tables")
dir_data <- file.path(dir_dst, "intermediate_data", plot_id)

dst_csv <- file.path(dir_csv, paste0(plot_id, "_ms_metrics.csv"))

# Stop if metrics have already been computed and recompute is not requested
if (file.exists(dst_csv) && !recompute) {
  cat("Metrics already computed for plot:", plot_id, "\n")
  quit(save = "no", status = 0)
}

# Load input data
if (!file.exists(src_bnd)) {
  stop(paste("The boundary file does not exist:", src_bnd))
}
if (!file.exists(src_treecover)) {
  stop(paste("The tree cover file does not exist:", src_treecover))
}

cat("Processing plot:", plot_id, "\n")

template <- terra::rast(
  file.path(dir_ms, "nir", plot_id, paste0(plot_id, "_nir.tif"))
  )
boundary <- terra::vect(src_bnd) %>%
  terra::project(template)
tree_cover <- terra::rast(src_treecover) %>%
  terra::project(template, method = "near") %>%
  terra::resample(template, method = "near") %>%
  terra::crop(boundary)

#>---------------------------------------------------------------------<
#> Functions
#>---------------------------------------------------------------------<
# Load index raster
get_index <- function(
    index, plot_id, dir_main = dir_ms, bound = boundary,
    tree_cover_raster = NULL, templ = template
    ) {
  print(paste("Computing index:", index, "for plot:", plot_id))
  f_rst <- file.path(
    dir_main, index, plot_id, paste0(plot_id, "_", index, ".tif")
    )
  if (!file.exists(f_rst)) {
    stop(paste("The index raster does not exist:", f_rst))
  }
  rst <- terra::rast(f_rst) %>%
    terra::crop(bound)
  
  if (!is.null(tree_cover_raster)) {
    rst <- terra::mask(
        rst, mask = tree_cover_raster, maskvalues = 0, updatevalue = udv
        )
  }

  return(terra::global(rst, fun = "mean", na.rm = TRUE)[1, 1])
}

get_ndvi2 <- function(
    plot_id, dir_main = dir_ms, bound = boundary, tree_cover_raster = NULL,
    templ = template
    ) {
  print(paste("Computing NDVI2 for plot:", plot_id))
  f_red <- file.path(
    dir_main, "red", plot_id, paste0(plot_id, "_red.tif")
    )
  f_nir <- file.path(
    dir_main, "nir", plot_id, paste0(plot_id, "_nir.tif")
    )
  if (!file.exists(f_red)) {
    stop(paste("The red raster does not exist:", f_red))
  }
  if (!file.exists(f_nir)) {
    stop(paste("The nir raster does not exist:", f_nir))
  }
  if (!is.null(tree_cover_raster)) {
    rst_red <- terra::rast(f_red) %>%
        terra::crop(bound) %>%
      terra::mask(
        mask = tree_cover_raster,
        maskvalues = 0, updatevalue = NA
        )
    rst_nir <- terra::rast(f_nir) %>%
      terra::crop(bound) %>%
      terra::mask(
        mask = tree_cover_raster,
        maskvalues = 0, updatevalue = NA
        )
  } else {
    rst_red <- terra::rast(f_red) %>%
      terra::crop(bound)
    rst_nir <- terra::rast(f_nir) %>%
      terra::crop(bound)
  }

  red <- terra::global(rst_red, fun = "mean", na.rm = TRUE)[1, 1]
  nir <- terra::global(rst_nir, fun = "mean", na.rm = TRUE)[1, 1]

  return((nir - red) / (nir + red))
}

#>---------------------------------------------------------------------<
#> Compute metrics
#>---------------------------------------------------------------------<
# Save metrics to CSV
df_metrics <- data.frame(
  plot_id = plot_id,
  lci_total = get_index("lci", plot_id, tree_cover_raster = NULL),
  lci_trees = get_index("lci", plot_id, tree_cover_raster = tree_cover),
  ndvi_total = get_index("ndvi", plot_id, tree_cover_raster = NULL),
  ndvi_trees = get_index("ndvi", plot_id, tree_cover_raster = tree_cover),
  ndvi2_total = get_ndvi2(plot_id, tree_cover_raster = NULL),
  ndvi2_trees = get_ndvi2(plot_id, tree_cover_raster = tree_cover),
  osavi_total = get_index("osavi", plot_id, tree_cover_raster = NULL),
  osavi_trees = get_index("osavi", plot_id, tree_cover_raster = tree_cover),
  gndvi_total = get_index("gndvi", plot_id, tree_cover_raster = NULL),
  gndvi_trees = get_index("gndvi", plot_id, tree_cover_raster = tree_cover),
  ndre_total = get_index("ndre", plot_id, tree_cover_raster = NULL),
  ndre_trees = get_index("ndre", plot_id, tree_cover_raster = tree_cover)
)

write.csv(df_metrics, dst_csv, row.names = FALSE)

cat("Metrics saved to:", dst_csv, "\n")
cat("\nFinished.")