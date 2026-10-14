library("dplyr")

dir_main <- "/storage/poppman/lai_metrics"
dir_metrics <- file.path(dir_main, "metric_tables")
dir_licor <- "/dme/eunomia/plotdata/measurments/ground/par/licor"

# Functions
get_lai <- function(plot_id, dir_licor) {
  f_lai <- file.path(dir_licor, plot_id, paste0(plot_id, ".csv"))
  if (!file.exists(f_lai)) {
    return(data.frame(plot_id = plot_id, lai = NA))
  }
  if (!"LAI" %in% names(read.csv(f_lai))) {
    warning(paste("The LAI column is missing in the file:", f_lai))
    return(data.frame(plot_id = plot_id, lai = NA))
  }
  lai_data <- read.csv(f_lai) %>%
    dplyr::pull("LAI") %>%
    mean(na.rm = TRUE) %>%
    data.frame(plot_id = plot_id, lai = .)
  return(lai_data)
}

# UAV-based metrics
files <- list.files(dir_metrics, pattern = ".csv$", full.names = TRUE)

files_lidar <- files[grepl("_lidar_", files)]
files_multispectral <- files[grepl("_ms_", files)]

lidar_table <- do.call(rbind, lapply(files_lidar, read.csv))
multispectral_table <- do.call(rbind, lapply(files_multispectral, read.csv))

licor_lai <- do.call(
  rbind, lapply(list.dirs(dir_licor, full.name = FALSE), get_lai, dir_licor = dir_licor)
  )

metrics_table <- dplyr::full_join(
  lidar_table, multispectral_table,
  by = "plot_id"
  )

write.csv(
  metrics_table,
  file.path(dir_main, "plotdata.csv"),
  row.names = FALSE
  )