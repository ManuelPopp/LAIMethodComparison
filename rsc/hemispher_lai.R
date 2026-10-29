if (!require(hemispheR)) {
  install.packages("hemispheR")
}

library("hemispheR")
library("dplyr")

dir_raw <- "/dme/eunomia/plotdata/measurments/ground/visible_spectrum/hemispheric/rawdata"
dir_selected <- "/dme/eunomia/plotdata/measurments/ground/visible_spectrum/hemispheric/hemisfer"
plot_folders <- list.dirs(dir_selected, recursive = FALSE, full.names = FALSE)

for (plot_folder in plot_folders) {
  plot_id <- basename(plot_folder)
  raw_files <- list.files(
    file.path(dir_raw, plot_id),
    pattern = "\\.JPG$", full.names = TRUE
    )
  if (length(raw_files) == 0) {
    warning(paste("No raw JPG files found for plot:", plot_id))
    next
  }

  selected_files <- list.files(
    file.path(dir_selected, plot_id),
    pattern = "\\.JPG$", full.names = FALSE
    )

  if (length(selected_files) == 0) {
    warning(paste("No selected JPG files found for plot:", plot_id))
    next
  }

  lais <- list()
  for (file in selected_files) {
    f_image <- file.path(raw_files, plot_id, file)
    img <- hemispheR::import_fisheye(
      f_image,
      channel = 3,
      circ.mask = list(xc = 4752, yc = 2672, rc = 5450),
      circular = FALSE,
      gamma = 2.2,
      stretch = FALSE,
      display = TRUE,
      message = TRUE
    )

    img.bw <- hemispheR::binarize_fisheye(
      img,
      method = "Otsu",
      zonal = FALSE,
      manual = NULL,
      display = TRUE,
      export = FALSE
    )

    gap.frac <- hemispheR::gapfrac_fisheye(
      img.bw,
      maxVZA = 90,
      lens = "FC-E8",
      startVZA = 0,
      endVZA = 70,
      nrings = 7,
      nseg = 8,
      display = TRUE,
      message = FALSE
    )

    lais[[file]] <- hemispheR::canopy_fisheye(gap.frac) %>%
      dplyr::rename(effective_LAI = Le, actual_LAI = L)
  }

  lai <- do.call(rbind, lais) %>%
    dplyr::mutate(plot_id = plot_id)
  write.csv(
    lai,
    file.path("/lud11/poppman/tmp/lai", plot_id, paste0(plot_id, "_hemisfer_lai.csv")),
    row.names = FALSE
    )
  )
}
