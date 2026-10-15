library("dplyr")
library("tidyr")
library("MASS")
library("ggplot2")

f_csv <- "L:/poppman/tmp/lai/plotdata.csv"

df <- read.csv(f_csv)

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
  dplyr::select(-lai)

results <- data.frame(
  variable = names(vars),
  correlation = sapply(vars, \(x) cor(x, df$lai, use = "complete.obs"))
) %>%
  dplyr::mutate(explained_variance = correlation^2) %>%
  dplyr::arrange(desc(explained_variance))

log_results <- data.frame(
  variable = names(vars),
  correlation = sapply(vars, \(x) {
    cor(log(x), df$lai, use = "complete.obs")
  })
) %>%
  dplyr::mutate(explained_variance = correlation^2) %>%
  dplyr::arrange(desc(explained_variance))

comparison <- results %>%
  dplyr::mutate(relationship = "linear") %>%
  dplyr::bind_rows(
    log_results %>%
      dplyr::mutate(relationship = "log-linear")
    ) %>%
  dplyr::arrange(desc(explained_variance))
rownames(comparison) <- NULL

comparison[which(grepl("trees", comparison$variable)),]

# Look at some correlations
plot(lai ~ log(chm_p90), data = df)
plot(lai ~ ndvi_trees, data = df)
plot(lai ~ log(lp1_median), data = df)
plot(lai ~ log(lp2_median), data = df)
plot(lai ~ log(abri_mean), data = df)

summary(lm(lai ~ tree_cover_fraction, data = df))
summary(lm(lai ~ log(chm_p90), data = df))
summary(lm(lai ~ ndvi_trees, data = df))
summary(lm(lai ~ ndvi_trees + chm_p90, data = df))

# Compare plots
gg_df <- df %>%
  dplyr::select(c("plot_id", "lai", comparison$variable[1:4])) %>%
  tidyr::pivot_longer(
    cols = -c(plot_id, lai),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    value = ifelse(
      metric %in% comparison$variable[comparison$relationship == "log-linear"],
      log(value), value
      ),
    metric = ifelse(
      metric %in% comparison$variable[comparison$relationship == "log-linear"],
      paste0(metric, " (log-transformed)"),
      metric
    ),
    metric = stringr::str_replace(
      metric, "^[^_]+", \(x) stringr::str_to_upper(x)
    ),
    metric = sub("_", " ", metric)
  )

ggplot2::ggplot(data = gg_df, ggplot2::aes(x = value, y = lai)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(method = "lm") +
  ggplot2::facet_wrap(. ~ metric, ncol = 2, scales = "free_x") +
  ggplot2::geom_text(
    data = gg_df %>%
      dplyr::group_by(metric) %>%
      dplyr::summarise(
        r2 = summary(lm(lai ~ value))$r.squared,
        value = Inf,
        lai = -Inf,
        .groups = "drop"
      ),
    ggplot2::aes(
      x = value, y = lai,
      label = paste0("R² = ", round(r2, 2))
    ),
    hjust = 1.1, vjust = -0.5,
    inherit.aes = FALSE
  ) +
  ggplot2::theme_classic() +
  ggplot2::theme(
    strip.background = ggplot2::element_blank()
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