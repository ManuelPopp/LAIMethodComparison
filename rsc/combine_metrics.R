dir_metrics <- "/storage/poppman/lai_metrics/metric_tables/"
files <- list.files(dir_metrics, pattern = ".csv$", full.names = TRUE)

files_lidar <- files[grepl("lidar", files)]
lidar_table <- do.call(rbind, lapply(files_lidar, read.csv))