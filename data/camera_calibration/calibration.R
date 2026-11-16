library("dplyr")

f_src <- "Calibration_measurements.csv"
f_dst <- "Calibration_coefficients.csv"

df <- read.csv(f_src, skip = 7, header = TRUE)

get_coefficients <- function(df, lens_name, r2_min = 0.99998) {
  r2 <- 0
  pol_deg <- 1
  while (r2 < r2_min & pol_deg <= 5) {
    mod <- subset(df, Lens == lens_name) %>%
      lm(RelativePixels ~ poly(
        Radians, degree = pol_deg, raw = TRUE) - 1,
        data = .
        )
    r2 <- summary(mod)$r.squared
    pol_deg <- pol_deg + 1
  }
  if (r2 < r2_min | pol_deg > 5) {
    stop("Failed fitting polynomial with minimum R2", r2_min)
  }
  coefs <- as.data.frame(t(coefficients(mod)))
  names(coefs) <- letters[1:ncol(coefs)]
  out_df <- cbind(
    data.frame(
      LensName = lens_name,
      Description = lens_name,
      MaxAngle = max(df$Degrees),
      Coefficients = 0
      ),
    coefs
    )
  return(out_df)
}

lapply(X = unique(df$Lens), FUN = get_coefficients, df = df) %>%
  dplyr::bind_rows() %>%
  write.csv(f_dst, row.names = FALSE)
