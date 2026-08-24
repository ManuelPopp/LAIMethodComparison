require("terra")
require("tidyterra")
require("dplyr")
require("readxl")
require("ggplot2")
library("foreach")
library("doParallel")

extract_values <- function(raster, locations, radius = 3) {
  circles <- terra::buffer(locations, width = radius)
  extracted <- terra::extract(
    raster, circles, fun = "mean", cells = FALSE, ID = FALSE, na.rm = TRUE
    )
  return(extracted)
}

cor.coef <- function(x, y) {
  coef <- cor.test(x, y, method = "pearson")$estimate
  return(coef)
}

p.value = function(x, y) {
  p <- cor.test(x, y, method = "pearson")$p.value
  return(p)
}

match_lai <- function(folder, index = "NDVI", radius = 10, mask_herbs = FALSE) {
  indices <- c(
    "LCI", "NDRE", "OSAVI", "NDVI", "GNDVI"
  )
  if (is.numeric(index)) {
    index <- indices[index]
  }
  f_rst <- file.path(folder, "DJITerra", paste0(index, ".tif"))
  f_lai <- file.path(
    folder, "Licor", paste0(basename(folder), "_Processed_coords.xlsx")
    )
  
  if (file.exists(f_rst) & file.exists(f_lai)) {
    print(paste0("File complete:", folder))
    rst <- terra::rast(f_rst) %>%
      terra::clamp(lower = 0, upper = 1, values = TRUE)
    
    if (mask_herbs) {
      dem <- terra::rast(
        file.path(folder, "DJITerra", paste0("dem", ".tif"))
        )
      dsm <- terra::rast(
        file.path(folder, "DJITerra", paste0("dsm", ".tif"))
      )
      mask <- (dsm - dem) >= 5
      rst <- terra::mask(rst, mask)
    }
    d_lai <- readxl::read_xlsx(
      path = f_lai
      )[, c("lon", "lat", "LAI", "LAI_File")] %>%
      dplyr::mutate(
        sample = tools::file_path_sans_ext(LAI_File)
        ) %>%
      dplyr::group_by(sample) %>%
      dplyr::summarise(
        lat = first(lat),
        lon = first(lon),
        lai = mean(LAI)
      )
    
    lai_vect <- terra::vect(
      d_lai, geom = c("lon", "lat"), crs = "epsg:4326"
      ) %>%
      terra::project(rst)
    
    index_vals <- extract_values(
      raster = rst,
      locations = lai_vect, radius = radius
    )
    
    lai_vals <- lai_vect$lai
    df <- data.frame(
      index_value = unname(index_vals), lai = lai_vals
    )
  } else {
    print(paste0("Files incomplete: ", folder))
    df <- data.frame(index_value = c(NA), lai = c(NA))
  }
  df$id <- basename(folder)
  return(df)
}

plots <- list.dirs("F:/FIELDWORK", recursive = FALSE)
metadata <- readxl::read_xlsx(
  path = "F:/FIELDWORK/plot_metadata.xlsx"
)

gridsearch_df <- data.frame()
indices <- c(
  "LCI", "NDRE", "OSAVI", "NDVI", "GNDVI"
)

cores <- parallel::detectCores()
cl <- parallel::makeCluster(cores - 1)
doParallel::registerDoParallel(cl)

gridsearch_df <- foreach::foreach(
  i = indices, .combine = "rbind", .packages = "dplyr"
  ) %:%
  foreach::foreach(
    r = seq(1, 26), .combine = "rbind", .packages = "dplyr"
    ) %dopar% {
    df <- do.call(
      rbind,
      lapply(X = plots, FUN = match_lai, index = i, radius = r)
    ) %>%
      dplyr::left_join(
        metadata, by = "id"
      )
    
    ct <- cor.test(
      df$lai,
      df$index_value,
      method = "pearson"
    )
    
    mod <- lm(lai ~ index_value, data = df)
    data.frame(
      radius = r,
      index = i,
      r_squared = summary(mod)$r.squared,
      p_value = ct$p.value
      )
}

parallel::stopCluster(cl)

ggplot2::ggplot(
  data = gridsearch_df,
  ggplot2::aes(x = radius, y = r_squared, colour = index)
) +
  ggplot2::geom_line() +
  ggplot2::geom_point() +
  ggplot2::theme_bw()

params <- gridsearch_df[
  which(gridsearch_df$r_squared == max(gridsearch_df$r_squared)),
  ]

df <- do.call(
  rbind,
  lapply(
    X = plots, FUN = match_lai, index = params$index, radius = params$radius
    )
) %>%
  dplyr::left_join(
    metadata, by = "id"
  )

# Measurement-level
ggplot2::ggplot(
  data = df, ggplot2::aes(x = lai, y = index_value, colour = vegetation_type)
) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(method = "lm") +
  ggplot2::geom_smooth(colour = "black", method = "lm") +
  ggplot2::theme_bw() +
  ggplot2::facet_wrap(.~vegetation_type)

# Stand-level
dfg <- df %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(
    index_value = median(index_value, na.rm = TRUE),
    lai = median(lai, na.rm = TRUE),
    vegetation_type = first(vegetation_type)
  )

ggplot2::ggplot(
  data = dfg,
  ggplot2::aes(x = index_value, y = lai, colour = vegetation_type)
) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(colour = "black", method = "lm") +
  ggplot2::theme_bw()

df_plot <- df %>%
  dplyr::group_by(id) %>%
  dplyr::summarise(
    index_value = mean(index_value, na.rm = TRUE),
    lai = mean(lai, na.rm = TRUE)
  )
