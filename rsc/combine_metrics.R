library("dplyr")

dir_main <- "/storage/poppman/lai_metrics"
dir_metrics <- file.path(dir_main, "metric_tables")
dir_licor <- "/dme/eunomia/plotdata/measurments/ground/par/licor"
dir_hemisfer <- "/dme/eunomia/plotdata/measurments/ground/visible_spectrum/hemispheric/lai"

# Functions
get_lai <- function(plot_id, dir_licor) {
  f_lai <- file.path(dir_licor, plot_id, paste0(plot_id, ".csv"))
  if (!file.exists(f_lai)) {
    return(data.frame(plot_id = tolower(plot_id), licor_lai = NA))
  }
  if (!"LAI" %in% names(read.csv(f_lai))) {
    warning(paste("The LAI column is missing in the file:", f_lai))
    return(data.frame(plot_id = tolower(plot_id), licor_lai = NA))
  }
  lai_data <- read.csv(f_lai) %>%
    dplyr::pull("LAI") %>%
    mean(na.rm = TRUE) %>%
    data.frame(plot_id = tolower(plot_id), licor_lai = .)
  return(lai_data)
}

get_hemlai <- function(plot_id, dir_hemisfer) {
  f_hlai <- file.path(dir_hemisfer, plot_id, paste0(plot_id, ".dat"))
  if (!file.exists(f_hlai)) {
    return(data.frame(plot_id = tolower(plot_id), hemisfer_lai = NA))
  }
  if (
    !"LAI2000" %in% names(
      read.table(
        f_hlai, sep = "\t", header = TRUE, fill = TRUE, quote = "",
        comment.char = ""
       )
      )
    ) {
    warning(paste("The LAI 2000 column is missing in the file:", f_hlai))
    return(data.frame(plot_id = tolower(plot_id), hemisfer_lai = NA))
  }
  lai_data <- base::suppressWarnings(
    read.table(
      f_hlai, sep = "\t", header = FALSE, fill = TRUE, quote = "",
      comment.char = ""
    )
   ) %>%
    dplyr::rename(Picture = V1, Correction = V3, LAI = V5) %>%
    dplyr::mutate(
      Picture = trimws(Picture),
      Correction = trimws(Correction),
      LAI = as.numeric(trimws(LAI))
      ) %>%
    dplyr::filter(
      Picture == "synthesis",
      Correction == "S. & al. (2007)"
      ) %>%
    dplyr::pull(LAI) %>%
    data.frame(plot_id = tolower(plot_id), hemisfer_lai = .)
  return(lai_data)
}

# UAV-based metrics
files <- list.files(dir_metrics, pattern = ".csv$", full.names = TRUE)

files_lidar <- files[grepl("_lidar_", files)]
files_multispectral <- files[grepl("_ms_", files)]

lidar_table <- do.call(rbind, lapply(files_lidar, read.csv))
multispectral_table <- do.call(rbind, lapply(files_multispectral, read.csv))

licor_lai <- do.call(
  rbind,
  lapply(list.dirs(dir_licor, full.name = FALSE), get_lai, dir_licor = dir_licor)
  )

hemisfer_lai <- do.call(
  rbind,
  lapply(
    list.dirs(dir_hemisfer, full.name = FALSE),
    get_hemlai, dir_hemisfer = dir_hemisfer
    )
  )

metrics_table <- dplyr::full_join(
  lidar_table, multispectral_table,
  by = "plot_id"
  ) %>%
    dplyr::full_join(licor_lai, by = "plot_id") %>%
    dplyr::full_join(hemisfer_lai, by = "plot_id")

write.csv(
  metrics_table,
  file.path(dir_main, "plotdata.csv"),
  row.names = FALSE
  )