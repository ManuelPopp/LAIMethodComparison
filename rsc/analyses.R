library("dplyr")
library("tidyr")
library("MASS")
library("ggplot2")

f_csv <- "L:/poppman/tmp/lai/plotdata.csv"
f_info <- "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/switchdrive/PhD/git/LAIMethodComparison/data/plot_info.csv"

df <- read.csv(f_csv) %>%
  dplyr::mutate(
    ndvi_cover = ndvi_trees * tree_cover_fraction
  )
df_info <- read.csv(f_info)

# Plot varibale correlations
df %>%
  dplyr::select(-plot_id) %>%
  stats::cor(use = "pairwise.complete.obs", method = "spearman") %>%
  corrplot::corrplot(
    method = "color", type = "upper",
    tl.col = "black", tl.cex = 0.7, cl.cex = 0.7,
    order = "hclust",
    col = rev(RColorBrewer::brewer.pal(n = 8, name = "RdBu"))
  )

# Get LAI correlations
vars <- df %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(-licor_lai)

calc_lm_metrics <- function(x, y, log_transform = FALSE) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  valid <- is.finite(x) & is.finite(y)
  
  if (log_transform) {
    valid <- valid & x > 0
  }
  
  x <- x[valid]
  y <- y[valid]
  
  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(c(R2adj = NA_real_, R2 = NA_real_))
  }
  
  if (log_transform) {
    x <- log(x)
  }
  
  tryCatch({
    model <- lm(y ~ x)
    c(
      R2adj = summary(model)$adj.r.squared,
      R2 = summary(model)$r.squared
    )
  }, error = function(e) {
    c(R2adj = NA_real_, R2 = NA_real_)
  })
}

results <- t(sapply(vars, \(x) {
  calc_lm_metrics(x, df$licor_lai)
})) |>
  as.data.frame() |>
  dplyr::mutate(variable = names(vars), .before = 1) |>
  dplyr::arrange(dplyr::desc(R2adj))

log_results <- t(sapply(vars, \(x) {
  calc_lm_metrics(x, df$licor_lai, log_transform = TRUE)
})) |>
  as.data.frame() |>
  dplyr::mutate(variable = names(vars), .before = 1) |>
  dplyr::arrange(dplyr::desc(R2adj))

comparison <- results %>%
  dplyr::mutate(relationship = "linear") %>%
  dplyr::bind_rows(
    log_results %>%
      dplyr::mutate(relationship = "log-linear")
    ) %>%
  dplyr::arrange(desc(R2adj))
rownames(comparison) <- NULL

# Look at some correlations
plot(licor_lai ~ hemisfer_lai, data = df)
plot(licor_lai ~ log(chm_p90), data = df)
plot(licor_lai ~ ndvi_trees, data = df)
plot(licor_lai ~ log(lp1_median), data = df)
plot(licor_lai ~ log(lp2_median), data = df)
plot(licor_lai ~ log(abri_mean), data = df)
plot(licor_lai ~ log(iri_mean), data = df)

summary(lm(licor_lai ~ hemisfer_lai, data = df))
summary(lm(licor_lai ~ tree_cover_fraction, data = df))
summary(lm(licor_lai ~ log(chm_p90), data = df))
summary(lm(licor_lai ~ ndvi_trees, data = df))
summary(lm(licor_lai ~ ndvi_trees + chm_p90, data = df))
summary(lm(licor_lai ~ ndvi_cover, data = df))

# Comparison plots
biomes <- df_info %>%
  dplyr::mutate(Plot = substr(PlotID, 1, 2)) %>%
  dplyr::group_by(Plot, Biome) %>%
  dplyr::summarise(
    plot = tolower(dplyr::first(Plot)),
    biome = vapply(
      strsplit(dplyr::first(Biome), " "),
      \(x) paste(head(x, 3), collapse = " "),
      character(1)
    ),
    .groups = "drop"
    )

top_vars <- comparison %>%
  dplyr::slice_head(n = 6) %>%
  dplyr::select(variable, relationship, R2, R2adj)

log_vars <- top_vars %>%
  dplyr::filter(relationship == "log-linear") %>%
  dplyr::pull(variable)

gg_df <- df %>%
  dplyr::select(c("plot_id", "licor_lai", top_vars$variable)) %>%
  tidyr::pivot_longer(
    cols = -c(plot_id, licor_lai),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    value = ifelse(metric %in% log_vars, log(value), value),
    metric = ifelse(
      metric %in% log_vars,
      paste0(metric, " (log-transformed)"),
      metric
    ),
    metric = stringr::str_replace(
      metric, "^[^_]+", \(x) stringr::str_to_upper(x)
    ),
    metric = sub("_", " ", metric),
    plot = substr(plot_id, 1, 2)
  ) %>%
  dplyr::left_join(y = biomes, by = "plot")

gg_mods <- ggplot2::ggplot(data = gg_df, ggplot2::aes(x = value, y = licor_lai)) +
  ggplot2::geom_point(ggplot2::aes(colour = biome, shape = biome)) +
  ggplot2::geom_smooth(method = "lm") +
  ggplot2::facet_wrap(. ~ metric, ncol = 3, scales = "free_x") +
  ggplot2::geom_text(
    data = gg_df %>%
      dplyr::group_by(metric) %>%
      dplyr::group_modify(\(d, key) {
        d <- d[is.finite(d$value) & is.finite(d$licor_lai), ]
        
        if (nrow(d) < 3 || dplyr::n_distinct(d$value) < 2) {
          return(data.frame(
            r2 = NA_real_,
            rmse = NA_real_,
            value = Inf,
            licor_lai = -Inf
          ))
        }
        
        model <- lm(licor_lai ~ value, data = d)
        
        data.frame(
          r2 = summary(model)$adj.r.squared,
          rmse = sqrt(mean(residuals(model)^2)),
          value = Inf,
          licor_lai = -Inf
        )
      }) %>%
      dplyr::ungroup(),
    ggplot2::aes(
      x = value, y = licor_lai,
      label = paste0(
        "R² = ", round(r2, 2),
        "\nRMSE = ", round(rmse, 2)
      )
    ),
    hjust = 1.1, vjust = -0.5,
    inherit.aes = FALSE
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    strip.background = ggplot2::element_blank(),
    legend.position = "bottom"
  ) +
  ggplot2::ylab("Leaf area index (Licor LAI 2200C)")

# Identify the "ideal model"----------------------------------------------------
# Identify numeric predictors
vars <- base::setdiff(names(df)[sapply(df, is.numeric)], "lai")

# Keep log terms only for strictly positive, finite variables
log_vars <- vars[sapply(vars, \(v) {
  x <- df[[v]]
  all(is.finite(x[!is.na(x)]) & x[!is.na(x)] > 0) &&
    any(!is.na(x))
})]

# Build candidate formula with linear and log terms
terms <- c(vars, paste0("log(", log_vars, ")"))

full_formula <- as.formula(
  paste("lai ~", paste(terms, collapse = " + "))
)

# Keep complete, finite observations for all candidate variables
model_vars <- c("lai", vars)

df_model <- df[complete.cases(df[, model_vars]), ]
df_model <- df_model[
  apply(df_model[, model_vars], 1, \(x) all(is.finite(x))),
]

# Fit null and full models on identical observations
null_model <- lm(lai ~ 1, data = df_model)
full_model <- lm(full_formula, data = df_model)

# Forward-backward stepwise selection using AIC
best_model <- stats::step(
  null_model,
  scope = list(
    lower = formula(null_model),
    upper = formula(full_model)
  ),
  direction = "both",
  trace = 0
)

# Results
summary(best_model)
formula(best_model)
AIC(best_model)

f_stillberg <- "D:/onedrive/OneDrive - Eidg. Forschungsanstalt WSL/shared/Stillberg/Summary_th1m.csv"
df_2 <- read.csv(f_stillberg) %>%
  dplyr::rename(
    ndvi_trees = tree_ndvi_mean,
    tree_cover_fraction = tree_cover,
    lai = LAI
  ) %>%
  dplyr::mutate(dataset = "Stillberg") %>%
  dplyr::select(
    ndvi_trees, tree_cover_fraction, lai, dataset
  ) %>%
  rbind(dplyr::select(df %>% dplyr::mutate(dataset = "Global"), ndvi_trees, tree_cover_fraction, lai, dataset)
  )

ggplot2::ggplot(
  data = df_2,
  ggplot2::aes(x = ndvi_trees * tree_cover_fraction, y = lai)
  ) +
  ggplot2::geom_point(ggplot2::aes(colour = dataset)) +
  ggplot2::geom_smooth(ggplot2::aes(colour = dataset), method = "lm") +
  ggplot2::geom_smooth(method = "lm")
