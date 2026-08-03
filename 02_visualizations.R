## Figures and tables for the manuscript

# load packages & functions
library(dplyr)
library(fmsb)
library(ggplot2)
library(openxlsx2)
library(patchwork)
library(scales)
library(scico)
library(tidyr)
library(xtable)

source("functions/problems.R")


# Help functions ---------------------------------------------------------------
scale_fill_diverging_fixed0 <- function(
    data = NULL,
    var = NULL,
    colors = c(
      "#001260", "#06558B", "#71A7C4",
      "#EBE5E0",
      "#D29773", "#AA4613", "#590007"
    ),
    name = NULL,
    limits = NULL,
    breaks = 5) {
  # Diverging fill scale with zero fixed at the midpoint.
  # data, var : data frame and variable used to infer limits if not provided
  # colors    : vector of colors (must be odd length; middle color maps to 0)
  # limits    : optional limits; defaults to symmetric range around 0

  stopifnot(length(colors) %% 2 == 1)

  # infer symmetric limits if not supplied
  if (is.null(limits)) {
    vals <- data[[deparse(substitute(var))]]
    max_abs <- max(abs(vals), na.rm = TRUE)
    limits <- c(-max_abs, max_abs)
  }

  mid <- (length(colors) + 1) / 2

  # anchor values with 0 exactly in the centre
  neg_vals <- seq(limits[1], 0, length.out = mid)
  pos_vals <- seq(0, limits[2], length.out = mid)

  scale_fill_gradientn(
    colors = colors,
    values = scales::rescale(c(neg_vals, pos_vals[-1]), from = limits),
    limits = limits,
    name = name,
    breaks = scales::breaks_pretty(n = breaks),
    oob = scales::squish
  )
}


scale_fill_positive_fixed0 <- function(
    data = NULL,
    var = NULL,
    colors = c("#EBE5E0", "#D29773", "#AA4613", "#590007"),
    name = NULL,
    limits = NULL,
    breaks = 5) {
  # Sequential fill scale for non-negative values with zero fixed at the lower
  # bound. Stops with an error if negative values are present.
  # data, var : data frame and variable used to infer limits if not provided
  # colors    : colors for increasing positive values

  # infer limits if not supplied
  if (is.null(limits)) {
    vals <- data[[deparse(substitute(var))]]
    limits <- range(vals, na.rm = TRUE)
  }

  if (limits[1] < 0) {
    stop("Negative values detected; use scale_fill_diverging_fixed0().")
  }

  scale_fill_gradientn(
    colors = colors,
    values = scales::rescale(c(0, limits[2]), from = limits),
    limits = limits,
    name = name,
    breaks = scales::breaks_pretty(n = breaks),
    oob = scales::squish
  )
}


# Preprocessing ----------------------------------------------------------------

## load results
res <- readRDS("results/res.rds")

## define labels
# algorithm
alg_labels <- c(
  alg_arima = "ARIMA",
  alg_moving_avg = "CMA",
  alg_flsa = "FLSA",
  alg_lowess = "LOWESS",
  alg_pelt = "PELT",
  alg_piecewise_reg = "PR",
  alg_gam = "TPRS"
)
# from scico palette "batlow"
# sample(scico::scico(7, palette = "batlow"))
alg_col <- c(
  "#D29243",
  "#818231",
  "#FDAC9C",
  "#134D61",
  "#001959",
  "#3B6C55",
  "#F9CCF9"
)

# distribution
dist_labels <- c(
  "norm" = "Normal",
  "lognorm" = "Log-normal"
)

# change
change_labels <- c(
  "No change", 
  "Change"
)

# estimands
est_labels <- c(
  "range" = "Range",
  "var" = "Variance",
  "madm" = "MADM",
  "n_cpts" = "NCP"
)
est_labels_short <- unname(est_labels)
est_labels_short[2] <- "Var"

# change pattern
pattern_labels <- c(
  "No change",
  "Jump",
  "Linear",
  "Quadratic",
  "Recurrent"
)

# sample size categories
nobs_col <- scico::scico(11, palette = "vikO")[5:1]

## preprocess results
res <- res %>%
  mutate(algorithm = factor(
    algorithm,
    levels = names(alg_labels),
    labels = unname(alg_labels)
  )) %>%
  mutate(change = if_else(
    grepl("nochange", problem), 0, 1)) %>%
  mutate(change = factor(
    change,
    levels = c(0,1),
    labels = change_labels
  )) %>%
  mutate(pattern = factor(
    pattern,
    levels = c("nochange", "jump", "linear", "quadratic", "recurrent"),
    labels = pattern_labels
  ))


# Figures ----------------------------------------------------------------------

#* Figure 1 --------------------------------------------------------------------
# Visualization of the five systematic change patterns over the data collection
# period.
nobs <- 1000
dmaxf <- 1
p1 <- 0.25
p2 <- 0.5
p3 <- 0.68
syserrdata <- data.frame(
  pattern = c(
    rep("No change", times = nobs),
    rep("Jump", times = nobs),
    rep("Linear", times = nobs),
    rep("Quadratic", times = nobs),
    rep("Recurrent", times = nobs)
  ),
  proportion = rep(seq(1, nobs, by = 1) / nobs, times = 5),
  syserr = c(
    f_simulate_norm_nochange(nobs = nobs, snr = 1)$syschange,
    f_simulate_norm_jump(nobs = nobs, snr = 1, dmaxf = dmaxf)$syschange,
    f_simulate_norm_linear(nobs = nobs, snr = 1, dmaxf = dmaxf)$syschange,
    f_simulate_norm_quadratic(nobs = nobs, snr = 1, dmaxf = dmaxf)$syschange,
    f_simulate_norm_recurrent(nobs = nobs, snr = 1, dmaxf = dmaxf)$syschange
  )
)
syserrdata <- syserrdata %>%
  mutate(pattern = factor(pattern, levels = c("No change", "Jump", "Linear", "Quadratic", "Recurrent")))
vlinedata <- data.frame(
  pattern = c(
    "Jump", "Linear", "Linear", "Quadratic", "Quadratic", "Quadratic",
    "Recurrent", "Recurrent", "Recurrent"
  ),
  vline_x = c(p1, p1, p2, p1, p2, p3, p1, p2, p3),
  annot_text = c("p1", "p1", "p2", "p1", "p2", "p3", "p1", "p2", "p3")
)
vlinedata <- vlinedata %>%
  mutate(pattern = factor(pattern, levels = c("No change", "Jump", "Linear", "Quadratic", "Recurrent")))
hlinedata <- data.frame(
  pattern = c("Jump", "Linear", "Quadratic", "Recurrent"),
  hline_y = c(1, 1, 1, 1),
  annot_text = c("dmax", "dmax", "dmax", "dmax")
)
hlinedata <- hlinedata %>%
  mutate(pattern = factor(pattern, levels = c("No change", "Jump", "Linear", "Quadratic", "Recurrent")))

ggplot(data = syserrdata, aes(x = proportion, y = syserr)) +
  geom_point() +
  theme_minimal() +
  labs(x = "Proportion of observations", y = "Systematic change") +
  facet_wrap(~pattern, ncol = 1) +
  geom_vline(data = vlinedata, aes(xintercept = vline_x), color = "#06558B", linetype = "dashed") +
  geom_hline(data = hlinedata, aes(yintercept = hline_y), color = "#AA4613", linetype = "dashed") +
  geom_text(data = vlinedata, aes(x = vline_x, y = 0.25, label = annot_text), color = "#06558B") +
  geom_text(data = hlinedata, aes(x = 0.05, y = hline_y, vjust = 1, label = annot_text), color = "#AA4613")
ggsave("results/figures/Fig1_syschange.pdf", width = 170, height = 121, units = "mm", bg = "white", dpi = 300)


#* Figure 2 --------------------------------------------------------------------
# Radar plots of scaled and inverted bias for quantifying systematic changes in
# normally distributed data.
layout_matrix <- matrix(c(
  rep(1, 4), # row 1: Sensitivity title
  2:5, # row 2: radar plots 1-4 Sensitivity
  6:9, # row 3: radar plots 5-7 Sensitivity
  rep(10, 4), # row 4: Specificity title
  11:14, # row 5: radar plots 1-4 Specificity
  15:18 # row 6: radar plots 5-7 Specificity
), nrow = 6, byrow = TRUE)
row_heights <- c(0.2, 1, 1, 0.2, 1, 1)

radar_norm_samp <- res %>%
  filter(distribution == "norm") %>%
  group_by(change, algorithm, nobs_cat) %>%
  summarise(
    Bias_range = mean(range_dif/sd_y, na.rm = TRUE),
    MSE_range = mean(range_dif^2/sd_y^2, na.rm = TRUE),
    Bias_variance = mean(var_dif/sd_y, na.rm = TRUE),
    MSE_variance = mean(var_dif^2/sd_y^2, na.rm = TRUE),
    Bias_MADM = mean(madm_dif/sd_y, na.rm = TRUE),
    MSE_MADM = mean(madm_dif^2/sd_y^2, na.rm = TRUE),
    Bias_NCP = mean(n_cpts_dif, na.rm = TRUE),
    MSE_NCP = mean(n_cpts_dif^2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # abs() is needed for bias and does not do harm for MSE
  mutate(
    across(
      matches("^(Bias|MSE)"),
      ~ 1 - scales::rescale(
        abs(.x),
        from = c(0, max(abs(.x), na.rm = TRUE))
      )
    )
  )

radar_norm_samp_bias <- radar_norm_samp %>%
  select(change, algorithm, nobs_cat, starts_with("Bias_"))

pdf("results/figures/Fig2_radar_norm_samp_bias.pdf",
  width = 6.69, #unit: inches (equals 170mm)
  height = 8 #maximum height: 225mm (8.85 in) for figure and legend
)
# Set up plotting area
layout(layout_matrix, heights = row_heights)
par(mar = c(1, 1, 1, 1), oma = c(1, 1, 1, 1), xpd = TRUE)

## Sensitivity
plot.new()
text(0.5, 0.5, "Change", cex = 1.4, font = 2)

# Loop through each method to plot
for (i in 1:length(alg_labels)) {
  alg <- unname(alg_labels[i])
  temp <- radar_norm_samp_bias %>%
    filter(algorithm == alg & change == "Change")
  radar_norm_alg <- rbind(rep(1, 4), 
                          rep(0, 4), 
                          temp[, c(-1, -2, -3)])
  rownames(radar_norm_alg) <- c("Max", "Min", as.character(temp$nobs_cat))

  fmsb::radarchart(
    radar_norm_alg,
    axistype = 1,
    vlabels = est_labels_short,
    pcol = nobs_col,
    plty = 1:5,
    plwd = 2,
    cglcol = "grey",
    cglty = 1,
    axislabcol = "black",
    caxislabels = c("0", "", "0.5", "", "1"),
    calcex = 0.9,
    vlcex = 0.9,
    title = "",
  )

  # Add method title
  mtext(unname(alg_labels[i]), side = 3, line = -0.3, cex = 0.8, font = 2)
}

## Specificity
plot.new() # to handle the missing 8. plot
plot.new()
text(0.5, 0.5, "No change", cex = 1.4, font = 2)

# Loop through each method to plot
for (i in 1:length(alg_labels)) {
  alg <- unname(alg_labels[i])
  temp <- radar_norm_samp_bias %>%
    filter(algorithm == alg & change == "No change")
  radar_norm_alg <- rbind(rep(1, 4), rep(0, 4), temp[, c(-1, -2, -3)])
  rownames(radar_norm_alg) <- c("Max", "Min", as.character(temp$nobs_cat))

  fmsb::radarchart(
    radar_norm_alg,
    axistype = 1,
    vlabels = est_labels_short,
    pcol = nobs_col,
    # pfcol = alpha(nobs_col, 0.2),
    plty = 1:5,
    plwd = 2,
    cglcol = "grey",
    cglty = 1,
    axislabcol = "black",
    caxislabels = c("0", "", "0.5", "", "1"),
    calcex = 0.9,
    vlcex = 0.9,
    title = "",
  )

  # Add method title
  mtext(unname(alg_labels[i]), side = 3, line = -0.3, cex = 0.8, font = 2)
}

# Add legend for sample sizes
plot.new()
legend("center",
  legend = unique(radar_norm_samp_bias$nobs_cat),
  col = nobs_col,
  lty = 1:5,
  lwd = 2,
  bty = "n",
  cex = 0.8,
  title = "Sample size",
  title.font = 2
)
dev.off()


#* Figure 3 --------------------------------------------------------------------
# Radar plots of scaled and inverted MSE for quantifying systematic changes in
# normally distributed data.
layout_matrix <- matrix(c(
  rep(1, 4), # row 1: Sensitivity title
  2:5, # row 2: radar plots 1-4 Sensitivity
  6:9, # row 3: radar plots 5-7 Sensitivity
  rep(10, 4), # row 4: Specificity title
  11:14, # row 5: radar plots 1-4 Specificity
  15:18 # row 6: radar plots 5-7 Specificity
), nrow = 6, byrow = TRUE)
row_heights <- c(0.2, 1, 1, 0.2, 1, 1)

radar_norm_samp <- res %>%
  filter(distribution == "norm") %>%
  group_by(change, algorithm, nobs_cat) %>%
  summarise(
    Bias_range = mean(range_dif/sd_y, na.rm = TRUE),
    MSE_range = mean(range_dif^2/sd_y^2, na.rm = TRUE),
    Bias_variance = mean(var_dif/sd_y, na.rm = TRUE),
    MSE_variance = mean(var_dif^2/sd_y^2, na.rm = TRUE),
    Bias_MADM = mean(madm_dif/sd_y, na.rm = TRUE),
    MSE_MADM = mean(madm_dif^2/sd_y^2, na.rm = TRUE),
    Bias_NCP = mean(n_cpts_dif, na.rm = TRUE),
    MSE_NCP = mean(n_cpts_dif^2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # abs() is needed for bias and does not do harm for MSE
  mutate(
    across(
      matches("^(Bias|MSE)"),
      ~ 1 - scales::rescale(
        abs(.x),
        from = c(0, max(abs(.x), na.rm = TRUE))
      )
    )
  )

radar_norm_samp_mse <- radar_norm_samp %>%
  select(change, algorithm, nobs_cat, starts_with("MSE_"))

pdf("results/figures/Fig3_radar_norm_samp_mse.pdf",
    width = 6.69, #unit: inches (equals 170mm)
    height = 8 #maximum height: 225mm (8.85 in) for figure and legend
)
# Set up plotting area
layout(layout_matrix, heights = row_heights)
par(mar = c(1, 1, 1, 1), oma = c(1, 1, 1, 1), xpd = TRUE)

## Sensitivity
plot.new()
text(0.5, 0.5, "Change", cex = 1.4, font = 2)

# Loop through each method to plot
for (i in 1:length(alg_labels)) {
  alg <- unname(alg_labels[i])
  temp <- radar_norm_samp_mse %>%
    filter(algorithm == alg & change == "Change")
  radar_norm_alg <- rbind(rep(1, 4), 
                          rep(0, 4), 
                          temp[, c(-1, -2, -3)])
  rownames(radar_norm_alg) <- c("Max", "Min", as.character(temp$nobs_cat))
  
  fmsb::radarchart(
    radar_norm_alg,
    axistype = 1,
    vlabels = est_labels_short,
    pcol = nobs_col,
    plty = 1:5,
    plwd = 2,
    cglcol = "grey",
    cglty = 1,
    axislabcol = "black",
    caxislabels = c("0", "", "0.5", "", "1"),
    calcex = 0.9,
    vlcex = 0.9,
    title = "",
  )
  
  # Add method title
  mtext(unname(alg_labels[i]), side = 3, line = -0.3, cex = 0.8, font = 2)
}

## Specificity
plot.new() # to handle the missing 8. plot
plot.new()
text(0.5, 0.5, "No change", cex = 1.4, font = 2)

# Loop through each method to plot
for (i in 1:length(alg_labels)) {
  alg <- unname(alg_labels[i])
  temp <- radar_norm_samp_mse %>%
    filter(algorithm == alg & change == "No change")
  radar_norm_alg <- rbind(rep(1, 4), rep(0, 4), temp[, c(-1, -2, -3)])
  rownames(radar_norm_alg) <- c("Max", "Min", as.character(temp$nobs_cat))
  
  fmsb::radarchart(
    radar_norm_alg,
    axistype = 1,
    vlabels = est_labels_short,
    pcol = nobs_col,
    # pfcol = alpha(nobs_col, 0.2),
    plty = 1:5,
    plwd = 2,
    cglcol = "grey",
    cglty = 1,
    axislabcol = "black",
    caxislabels = c("0", "", "0.5", "", "1"),
    calcex = 0.9,
    vlcex = 0.9,
    title = "",
  )
  
  # Add method title
  mtext(unname(alg_labels[i]), side = 3, line = -0.3, cex = 0.8, font = 2)
}

# Add legend for sample sizes
plot.new()
legend("center",
       legend = unique(radar_norm_samp_mse$nobs_cat),
       col = nobs_col,
       lty = 1:5,
       lwd = 2,
       bty = "n",
       cex = 0.8,
       title = "Sample size",
       title.font = 2
)
dev.off()


# Tables -----------------------------------------------------------------------

#* Table 2 & S1 ----------------------------------------------------------------
# Sensitivity: Bias and mean squared error for quantifying systematic changes in
# normally distributed data.

# by sample size category
sum_norm_bias_mse <- res %>%
  filter(distribution == "norm" & change == "Change") %>%
  group_by(nobs_cat, algorithm) %>%
  summarise(
    n = n(),
    Bias_range = mean(range_dif/sd_y),
    Bias_range_se = sd(range_dif/sd_y)/sqrt(n()),
    MSE_range = mean(range_dif^2/sd_y^2),
    MSE_range_se = sd(range_dif^2/sd_y^2)/sqrt(n()),
    Bias_variance = mean(var_dif/sd_y),
    Bias_variance_se = sd(var_dif/sd_y)/sqrt(n()),
    MSE_variance = mean(var_dif^2/sd_y^2),
    MSE_variance_se = sd(var_dif^2/sd_y^2)/sqrt(n()),
    Bias_MADM = mean(madm_dif/sd_y),
    Bias_MADM_se = sd(madm_dif/sd_y)/sqrt(n()),
    MSE_MADM = mean(madm_dif^2/sd_y^2),
    MSE_MADM_se = sd(madm_dif^2/sd_y^2)/sqrt(n()),
    Bias_NCP = mean(n_cpts_dif),
    Bias_NCP_se = sd(n_cpts_dif)/sqrt(n()),
    MSE_NCP = mean(n_cpts_dif^2),
    MSE_NCP_se = sd(n_cpts_dif^2)/sqrt(n()),
    .groups = "drop"
  )

# aggregated across sample size categories
sum_norm_bias_mse_agg <- res %>%
  filter(distribution == "norm" & change == "Change") %>%
  group_by(algorithm) %>%
  summarise(
    n = n(),
    Bias_range = mean(range_dif/sd_y),
    Bias_range_se = sd(range_dif/sd_y)/sqrt(n()),
    MSE_range = mean(range_dif^2/sd_y^2),
    MSE_range_se = sd(range_dif^2/sd_y^2)/sqrt(n()),
    Bias_variance = mean(var_dif/sd_y),
    Bias_variance_se = sd(var_dif/sd_y)/sqrt(n()),
    MSE_variance = mean(var_dif^2/sd_y^2),
    MSE_variance_se = sd(var_dif^2/sd_y^2)/sqrt(n()),
    Bias_MADM = mean(madm_dif/sd_y),
    Bias_MADM_se = sd(madm_dif/sd_y)/sqrt(n()),
    MSE_MADM = mean(madm_dif^2/sd_y^2),
    MSE_MADM_se = sd(madm_dif^2/sd_y^2)/sqrt(n()),
    Bias_NCP = mean(n_cpts_dif),
    Bias_NCP_se = sd(n_cpts_dif)/sqrt(n()),
    MSE_NCP = mean(n_cpts_dif^2),
    MSE_NCP_se = sd(n_cpts_dif^2)/sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(nobs_cat = "30-1000")

# combine
sum_norm_bias_mse <- sum_norm_bias_mse %>%
  mutate(
    nobs_cat = forcats::fct_expand(nobs_cat, "30–1000")
  ) %>%
  bind_rows(sum_norm_bias_mse_agg)

# compute borda ranking as a summary measure
sum_norm_bias_mse <- sum_norm_bias_mse %>%  
  group_by(nobs_cat) %>%
  mutate(across(
    matches("^(Bias|MSE)") & !ends_with("_se"),
    ~ rank(abs(.x), ties.method = "average"),
    # abs() is needed for bias and does not do harm for MSE
    .names = "{.col}_rank"
  )) %>%
  mutate(Bias_borda = rowSums(across(matches("^Bias.*_rank$"))),
         Bias_borda_scale = (Bias_borda - 4*1) / (4*7 - 4*1), 
         MSE_borda = rowSums(across(matches("^MSE.*_rank$"))),
         MSE_borda_scale = (MSE_borda - 4*1) / (4*7 - 4*1) ) %>% 
  ungroup() %>%
  select(-c(ends_with("_rank"), "Bias_borda", "MSE_borda"))

# highlight best values for latex output in bold
sum_norm_bias_mse_highlighted <- sum_norm_bias_mse %>%
  rename(
    Sample_size = nobs_cat,
    Method = algorithm
  ) %>%
  group_by(Sample_size) %>%
  mutate(across(
    .cols = matches("^(Bias|MSE)") & !ends_with("_se"),
    .fns = ~ {
      orig_vals <- round(.x, 2) # for display
      abs_vals <- abs(orig_vals) # for highlighting
      
      # smaller-is-better columns
      target_val <- min(abs_vals, na.rm = TRUE)
      is_target <- abs_vals == target_val

      # highlight original values (not absolute!)
      ifelse(
        is_target,
        paste0("\\textbf{", sprintf("%.2f", orig_vals), "}"),
        sprintf("%.2f", orig_vals)
      )
    }
  )) %>%
  ungroup() %>%
  mutate(
    Sample_size = as.character(Sample_size),
    Sample_size = replace(Sample_size, duplicated(Sample_size), "")
  )

## Table 2 (Bias)
sum_norm_bias_highlighted <- sum_norm_bias_mse_highlighted %>% 
  select(-starts_with("MSE"))

latex_table2_norm <- xtable::xtable(sum_norm_bias_highlighted, digits = 2)
print(latex_table2_norm,
  type = "latex",
  include.rownames = FALSE,
  booktabs = TRUE,
  sanitize.text.function = identity,
  file = "results/tables/table2_norm_bias_sens.tex"
)

## Table S1 (MSE)
sum_norm_mse_highlighted <- sum_norm_bias_mse_highlighted %>% 
  select(-starts_with("Bias"))

latex_tableS1_norm <- xtable::xtable(sum_norm_mse_highlighted, digits = 2)
print(latex_tableS1_norm,
      type = "latex",
      include.rownames = FALSE,
      booktabs = TRUE,
      sanitize.text.function = identity,
      file = "results/tables/tableS1_norm_mse_sens.tex"
)


#* Table 3 & S2 ----------------------------------------------------------------
# Specificity: Bias and mean squared error for quantifying systematic changes in 
# normally distributed data.

# by sample size category
sum_norm_bias_mse <- res %>%
  filter(distribution == "norm" & change == "No change") %>%
  group_by(nobs_cat, algorithm) %>%
  summarise(
    n = n(),
    Bias_range = mean(range_dif/sd_y),
    Bias_range_se = sd(range_dif/sd_y)/sqrt(n()),
    MSE_range = mean(range_dif^2/sd_y^2),
    MSE_range_se = sd(range_dif^2/sd_y^2)/sqrt(n()),
    Bias_variance = mean(var_dif/sd_y),
    Bias_variance_se = sd(var_dif/sd_y)/sqrt(n()),
    MSE_variance = mean(var_dif^2/sd_y^2),
    MSE_variance_se = sd(var_dif^2/sd_y^2)/sqrt(n()),
    Bias_MADM = mean(madm_dif/sd_y),
    Bias_MADM_se = sd(madm_dif/sd_y)/sqrt(n()),
    MSE_MADM = mean(madm_dif^2/sd_y^2),
    MSE_MADM_se = sd(madm_dif^2/sd_y^2)/sqrt(n()),
    Bias_NCP = mean(n_cpts_dif),
    Bias_NCP_se = sd(n_cpts_dif)/sqrt(n()),
    MSE_NCP = mean(n_cpts_dif^2),
    MSE_NCP_se = sd(n_cpts_dif^2)/sqrt(n()),
    .groups = "drop"
  )

# aggregated across sample size categories
sum_norm_bias_mse_agg <- res %>%
  filter(distribution == "norm" & change == "No change") %>%
  group_by(algorithm) %>%
  summarise(
    n = n(),
    Bias_range = mean(range_dif/sd_y),
    Bias_range_se = sd(range_dif/sd_y)/sqrt(n()),
    MSE_range = mean(range_dif^2/sd_y^2),
    MSE_range_se = sd(range_dif^2/sd_y^2)/sqrt(n()),
    Bias_variance = mean(var_dif/sd_y),
    Bias_variance_se = sd(var_dif/sd_y)/sqrt(n()),
    MSE_variance = mean(var_dif^2/sd_y^2),
    MSE_variance_se = sd(var_dif^2/sd_y^2)/sqrt(n()),
    Bias_MADM = mean(madm_dif/sd_y),
    Bias_MADM_se = sd(madm_dif/sd_y)/sqrt(n()),
    MSE_MADM = mean(madm_dif^2/sd_y^2),
    MSE_MADM_se = sd(madm_dif^2/sd_y^2)/sqrt(n()),
    Bias_NCP = mean(n_cpts_dif),
    Bias_NCP_se = sd(n_cpts_dif)/sqrt(n()),
    MSE_NCP = mean(n_cpts_dif^2),
    MSE_NCP_se = sd(n_cpts_dif^2)/sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(nobs_cat = "30-1000")

# combine
sum_norm_bias_mse <- sum_norm_bias_mse %>%
  mutate(
    nobs_cat = forcats::fct_expand(nobs_cat, "30–1000")
  ) %>%
  bind_rows(sum_norm_bias_mse_agg)

# compute borda ranking as a summary measure
sum_norm_bias_mse <- sum_norm_bias_mse %>%  
  group_by(nobs_cat) %>%
  mutate(across(
    matches("^(Bias|MSE)") & !ends_with("_se"),
    ~ rank(abs(.x), ties.method = "average"),
    # abs() is needed for bias and does not do harm for MSE
    .names = "{.col}_rank"
  )) %>%
  mutate(Bias_borda = rowSums(across(matches("^Bias.*_rank$"))),
         Bias_borda_scale = (Bias_borda - 4*1) / (4*7 - 4*1), 
         MSE_borda = rowSums(across(matches("^MSE.*_rank$"))),
         MSE_borda_scale = (MSE_borda - 4*1) / (4*7 - 4*1) ) %>% 
  ungroup() %>%
  select(-c(ends_with("_rank"), "Bias_borda", "MSE_borda"))

# highlight best values for latex output in bold
sum_norm_bias_mse_highlighted <- sum_norm_bias_mse %>%
  rename(
    Sample_size = nobs_cat,
    Method = algorithm
  ) %>%
  group_by(Sample_size) %>%
  mutate(across(
    .cols = matches("^(Bias|MSE)") & !ends_with("_se"),
    .fns = ~ {
      orig_vals <- round(.x, 2) # for display
      abs_vals <- abs(orig_vals) # for highlighting
      
      # smaller-is-better columns
      target_val <- min(abs_vals, na.rm = TRUE)
      is_target <- abs_vals == target_val
      
      # highlight original values (not absolute!)
      ifelse(
        is_target,
        paste0("\\textbf{", sprintf("%.2f", orig_vals), "}"),
        sprintf("%.2f", orig_vals)
      )
    }
  )) %>%
  ungroup() %>%
  mutate(
    Sample_size = as.character(Sample_size),
    Sample_size = replace(Sample_size, duplicated(Sample_size), "")
  )

## Table 3 (Bias)
sum_norm_bias_highlighted <- sum_norm_bias_mse_highlighted %>% 
  select(-starts_with("MSE"))

latex_table3_norm <- xtable::xtable(sum_norm_bias_highlighted, digits = 2)
print(latex_table3_norm,
      type = "latex",
      include.rownames = FALSE,
      booktabs = TRUE,
      sanitize.text.function = identity,
      file = "results/tables/table3_norm_bias_spec.tex"
)

## Table S2 (MSE)
sum_norm_mse_highlighted <- sum_norm_bias_mse_highlighted %>% 
  select(-starts_with("Bias"))

latex_tableS2_norm <- xtable::xtable(sum_norm_mse_highlighted, digits = 2)
print(latex_tableS2_norm,
      type = "latex",
      include.rownames = FALSE,
      booktabs = TRUE,
      sanitize.text.function = identity,
      file = "results/tables/tableS2_norm_mse_spec.tex"
)


# Suppl. tables ----------------------------------------------------------------

#* Table  S3 & S5 --------------------------------------------------------------
# Sensitivity: Bias and mean squared error for quantifying systematic changes in
# log-normally distributed data.

# by sample size category
sum_lognorm_bias_mse <- res %>%
  filter(distribution == "lognorm" & change == "Change") %>%
  group_by(nobs_cat, algorithm) %>%
  summarise(
    n = n(),
    Bias_range = mean(range_dif/sd_y),
    Bias_range_se = sd(range_dif/sd_y)/sqrt(n()),
    MSE_range = mean(range_dif^2/sd_y^2),
    MSE_range_se = sd(range_dif^2/sd_y^2)/sqrt(n()),
    Bias_variance = mean(var_dif/sd_y),
    Bias_variance_se = sd(var_dif/sd_y)/sqrt(n()),
    MSE_variance = mean(var_dif^2/sd_y^2),
    MSE_variance_se = sd(var_dif^2/sd_y^2)/sqrt(n()),
    Bias_MADM = mean(madm_dif/sd_y),
    Bias_MADM_se = sd(madm_dif/sd_y)/sqrt(n()),
    MSE_MADM = mean(madm_dif^2/sd_y^2),
    MSE_MADM_se = sd(madm_dif^2/sd_y^2)/sqrt(n()),
    Bias_NCP = mean(n_cpts_dif),
    Bias_NCP_se = sd(n_cpts_dif)/sqrt(n()),
    MSE_NCP = mean(n_cpts_dif^2),
    MSE_NCP_se = sd(n_cpts_dif^2)/sqrt(n()),
    .groups = "drop"
  )

# aggregated across sample size categories
sum_lognorm_bias_mse_agg <- res %>%
  filter(distribution == "lognorm" & change == "Change") %>%
  group_by(algorithm) %>%
  summarise(
    n = n(),
    Bias_range = mean(range_dif/sd_y),
    Bias_range_se = sd(range_dif/sd_y)/sqrt(n()),
    MSE_range = mean(range_dif^2/sd_y^2),
    MSE_range_se = sd(range_dif^2/sd_y^2)/sqrt(n()),
    Bias_variance = mean(var_dif/sd_y),
    Bias_variance_se = sd(var_dif/sd_y)/sqrt(n()),
    MSE_variance = mean(var_dif^2/sd_y^2),
    MSE_variance_se = sd(var_dif^2/sd_y^2)/sqrt(n()),
    Bias_MADM = mean(madm_dif/sd_y),
    Bias_MADM_se = sd(madm_dif/sd_y)/sqrt(n()),
    MSE_MADM = mean(madm_dif^2/sd_y^2),
    MSE_MADM_se = sd(madm_dif^2/sd_y^2)/sqrt(n()),
    Bias_NCP = mean(n_cpts_dif),
    Bias_NCP_se = sd(n_cpts_dif)/sqrt(n()),
    MSE_NCP = mean(n_cpts_dif^2),
    MSE_NCP_se = sd(n_cpts_dif^2)/sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(nobs_cat = "30-1000")

# combine
sum_lognorm_bias_mse <- sum_lognorm_bias_mse %>%
  mutate(
    nobs_cat = forcats::fct_expand(nobs_cat, "30–1000")
  ) %>%
  bind_rows(sum_norm_bias_mse_agg)

# compute borda ranking as a summary measure
sum_lognorm_bias_mse <- sum_lognorm_bias_mse %>%  
  group_by(nobs_cat) %>%
  mutate(across(
    matches("^(Bias|MSE)") & !ends_with("_se"),
    ~ rank(abs(.x), ties.method = "average"),
    # abs() is needed for bias and does not do harm for MSE
    .names = "{.col}_rank"
  )) %>%
  mutate(Bias_borda = rowSums(across(matches("^Bias.*_rank$"))),
         Bias_borda_scale = (Bias_borda - 4*1) / (4*7 - 4*1), 
         MSE_borda = rowSums(across(matches("^MSE.*_rank$"))),
         MSE_borda_scale = (MSE_borda - 4*1) / (4*7 - 4*1) ) %>% 
  ungroup() %>%
  select(-c(ends_with("_rank"), "Bias_borda", "MSE_borda"))

# highlight best values for latex output in bold
sum_lognorm_bias_mse_highlighted <- sum_lognorm_bias_mse %>%
  rename(
    Sample_size = nobs_cat,
    Method = algorithm
  ) %>%
  group_by(Sample_size) %>%
  mutate(across(
    .cols = matches("^(Bias|MSE)") & !ends_with("_se"),
    .fns = ~ {
      orig_vals <- round(.x, 2) # for display
      abs_vals <- abs(orig_vals) # for highlighting
      
      # smaller-is-better columns
      target_val <- min(abs_vals, na.rm = TRUE)
      is_target <- abs_vals == target_val
      
      # highlight original values (not absolute!)
      ifelse(
        is_target,
        paste0("\\textbf{", sprintf("%.2f", orig_vals), "}"),
        sprintf("%.2f", orig_vals)
      )
    }
  )) %>%
  ungroup() %>%
  mutate(
    Sample_size = as.character(Sample_size),
    Sample_size = replace(Sample_size, duplicated(Sample_size), "")
  )

## Table S3 (Bias)
sum_lognorm_bias_highlighted <- sum_lognorm_bias_mse_highlighted %>% 
  select(-starts_with("MSE"))

latex_tableS3_lognorm <- xtable::xtable(sum_lognorm_bias_highlighted, digits = 2)
print(latex_tableS3_lognorm,
      type = "latex",
      include.rownames = FALSE,
      booktabs = TRUE,
      sanitize.text.function = identity,
      file = "results/tables/tableS3_lognorm_bias_sens.tex"
)

## Table S5 (MSE)
sum_lognorm_mse_highlighted <- sum_lognorm_bias_mse_highlighted %>% 
  select(-starts_with("Bias"))

latex_tableS5_lognorm <- xtable::xtable(sum_lognorm_mse_highlighted, digits = 2)
print(latex_tableS5_lognorm,
      type = "latex",
      include.rownames = FALSE,
      booktabs = TRUE,
      sanitize.text.function = identity,
      file = "results/tables/tableS5_lognorm_mse_sens.tex"
)


#* Table S4 & S6 ---------------------------------------------------------------
# Specificity: Bias and mean squared error for quantifying systematic changes in
# log-normally distributed data.

# by sample size category
sum_lognorm_bias_mse <- res %>%
  filter(distribution == "lognorm" & change == "No change") %>%
  group_by(nobs_cat, algorithm) %>%
  summarise(
    n = n(),
    Bias_range = mean(range_dif/sd_y),
    Bias_range_se = sd(range_dif/sd_y)/sqrt(n()),
    MSE_range = mean(range_dif^2/sd_y^2),
    MSE_range_se = sd(range_dif^2/sd_y^2)/sqrt(n()),
    Bias_variance = mean(var_dif/sd_y),
    Bias_variance_se = sd(var_dif/sd_y)/sqrt(n()),
    MSE_variance = mean(var_dif^2/sd_y^2),
    MSE_variance_se = sd(var_dif^2/sd_y^2)/sqrt(n()),
    Bias_MADM = mean(madm_dif/sd_y),
    Bias_MADM_se = sd(madm_dif/sd_y)/sqrt(n()),
    MSE_MADM = mean(madm_dif^2/sd_y^2),
    MSE_MADM_se = sd(madm_dif^2/sd_y^2)/sqrt(n()),
    Bias_NCP = mean(n_cpts_dif),
    Bias_NCP_se = sd(n_cpts_dif)/sqrt(n()),
    MSE_NCP = mean(n_cpts_dif^2),
    MSE_NCP_se = sd(n_cpts_dif^2)/sqrt(n()),
    .groups = "drop"
  )

# aggregated across sample size categories
sum_lognorm_bias_mse_agg <- res %>%
  filter(distribution == "lognorm" & change == "No change") %>%
  group_by(algorithm) %>%
  summarise(
    n = n(),
    Bias_range = mean(range_dif/sd_y),
    Bias_range_se = sd(range_dif/sd_y)/sqrt(n()),
    MSE_range = mean(range_dif^2/sd_y^2),
    MSE_range_se = sd(range_dif^2/sd_y^2)/sqrt(n()),
    Bias_variance = mean(var_dif/sd_y),
    Bias_variance_se = sd(var_dif/sd_y)/sqrt(n()),
    MSE_variance = mean(var_dif^2/sd_y^2),
    MSE_variance_se = sd(var_dif^2/sd_y^2)/sqrt(n()),
    Bias_MADM = mean(madm_dif/sd_y),
    Bias_MADM_se = sd(madm_dif/sd_y)/sqrt(n()),
    MSE_MADM = mean(madm_dif^2/sd_y^2),
    MSE_MADM_se = sd(madm_dif^2/sd_y^2)/sqrt(n()),
    Bias_NCP = mean(n_cpts_dif),
    Bias_NCP_se = sd(n_cpts_dif)/sqrt(n()),
    MSE_NCP = mean(n_cpts_dif^2),
    MSE_NCP_se = sd(n_cpts_dif^2)/sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(nobs_cat = "30-1000")

# combine
sum_lognorm_bias_mse <- sum_lognorm_bias_mse %>%
  mutate(
    nobs_cat = forcats::fct_expand(nobs_cat, "30–1000")
  ) %>%
  bind_rows(sum_norm_bias_mse_agg)

# compute borda ranking as a summary measure
sum_lognorm_bias_mse <- sum_lognorm_bias_mse %>%  
  group_by(nobs_cat) %>%
  mutate(across(
    matches("^(Bias|MSE)") & !ends_with("_se"),
    ~ rank(abs(.x), ties.method = "average"),
    # abs() is needed for bias and does not do harm for MSE
    .names = "{.col}_rank"
  )) %>%
  mutate(Bias_borda = rowSums(across(matches("^Bias.*_rank$"))),
         Bias_borda_scale = (Bias_borda - 4*1) / (4*7 - 4*1), 
         MSE_borda = rowSums(across(matches("^MSE.*_rank$"))),
         MSE_borda_scale = (MSE_borda - 4*1) / (4*7 - 4*1) ) %>% 
  ungroup() %>%
  select(-c(ends_with("_rank"), "Bias_borda", "MSE_borda"))

# highlight best values for latex output in bold
sum_lognorm_bias_mse_highlighted <- sum_lognorm_bias_mse %>%
  rename(
    Sample_size = nobs_cat,
    Method = algorithm
  ) %>%
  group_by(Sample_size) %>%
  mutate(across(
    .cols = matches("^(Bias|MSE)") & !ends_with("_se"),
    .fns = ~ {
      orig_vals <- round(.x, 2) # for display
      abs_vals <- abs(orig_vals) # for highlighting
      
      # smaller-is-better columns
      target_val <- min(abs_vals, na.rm = TRUE)
      is_target <- abs_vals == target_val
      
      # highlight original values (not absolute!)
      ifelse(
        is_target,
        paste0("\\textbf{", sprintf("%.2f", orig_vals), "}"),
        sprintf("%.2f", orig_vals)
      )
    }
  )) %>%
  ungroup() %>%
  mutate(
    Sample_size = as.character(Sample_size),
    Sample_size = replace(Sample_size, duplicated(Sample_size), "")
  )

## Table S4 (Bias)
sum_lognorm_bias_highlighted <- sum_lognorm_bias_mse_highlighted %>% 
  select(-starts_with("MSE"))

latex_tableS4_lognorm <- xtable::xtable(sum_lognorm_bias_highlighted, digits = 2)
print(latex_tableS4_lognorm,
      type = "latex",
      include.rownames = FALSE,
      booktabs = TRUE,
      sanitize.text.function = identity,
      file = "results/tables/tableS4_lognorm_bias_spec.tex"
)

## Table S6 (MSE)
sum_lognorm_mse_highlighted <- sum_lognorm_bias_mse_highlighted %>% 
  select(-starts_with("Bias"))

latex_tableS6_norm <- xtable::xtable(sum_norm_mse_highlighted, digits = 2)
print(latex_tableS6_norm,
      type = "latex",
      include.rownames = FALSE,
      booktabs = TRUE,
      sanitize.text.function = identity,
      file = "results/tables/tableS6_norm_mse_spec.tex"
)


# Suppl. figures ---------------------------------------------------------------
#* Figure S1 --------------------------------------------------------------------
# Radar plots of scaled and inverted bias for quantifying systematic changes in
# log-normally distributed data.
layout_matrix <- matrix(c(
  rep(1, 4), # row 1: Sensitivity title
  2:5, # row 2: radar plots 1-4 Sensitivity
  6:9, # row 3: radar plots 5-7 Sensitivity
  rep(10, 4), # row 4: Specificity title
  11:14, # row 5: radar plots 1-4 Specificity
  15:18 # row 6: radar plots 5-7 Specificity
), nrow = 6, byrow = TRUE)
row_heights <- c(0.2, 1, 1, 0.2, 1, 1)

radar_lognorm_samp <- res %>%
  filter(distribution == "lognorm") %>%
  group_by(change, algorithm, nobs_cat) %>%
  summarise(
    Bias_range = mean(range_dif/sd_y, na.rm = TRUE),
    MSE_range = mean(range_dif^2/sd_y^2, na.rm = TRUE),
    Bias_variance = mean(var_dif/sd_y, na.rm = TRUE),
    MSE_variance = mean(var_dif^2/sd_y^2, na.rm = TRUE),
    Bias_MADM = mean(madm_dif/sd_y, na.rm = TRUE),
    MSE_MADM = mean(madm_dif^2/sd_y^2, na.rm = TRUE),
    Bias_NCP = mean(n_cpts_dif, na.rm = TRUE),
    MSE_NCP = mean(n_cpts_dif^2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # abs() is needed for bias and does not do harm for MSE
  mutate(
    across(
      matches("^(Bias|MSE)"),
      ~ 1 - scales::rescale(
        abs(.x),
        from = c(0, max(abs(.x), na.rm = TRUE))
      )
    )
  )

radar_lognorm_samp_bias <- radar_lognorm_samp %>%
  select(change, algorithm, nobs_cat, starts_with("Bias_"))

pdf("results/figures/FigS1_radar_lognorm_samp_bias.pdf",
    width = 6.69, #unit: inches (equals 170mm)
    height = 8 #maximum height: 225mm (8.85 in) for figure and legend
)
# Set up plotting area
layout(layout_matrix, heights = row_heights)
par(mar = c(1, 1, 1, 1), oma = c(1, 1, 1, 1), xpd = TRUE)

## Sensitivity
plot.new()
text(0.5, 0.5, "Change", cex = 1.4, font = 2)

# Loop through each method to plot
for (i in 1:length(alg_labels)) {
  alg <- unname(alg_labels[i])
  temp <- radar_lognorm_samp_bias %>%
    filter(algorithm == alg & change == "Change")
  radar_lognorm_alg <- rbind(rep(1, 4), 
                          rep(0, 4), 
                          temp[, c(-1, -2, -3)])
  rownames(radar_lognorm_alg) <- c("Max", "Min", as.character(temp$nobs_cat))
  
  fmsb::radarchart(
    radar_lognorm_alg,
    axistype = 1,
    vlabels = est_labels_short,
    pcol = nobs_col,
    plty = 1:5,
    plwd = 2,
    cglcol = "grey",
    cglty = 1,
    axislabcol = "black",
    caxislabels = c("0", "", "0.5", "", "1"),
    calcex = 0.9,
    vlcex = 0.9,
    title = "",
  )
  
  # Add method title
  mtext(unname(alg_labels[i]), side = 3, line = -0.3, cex = 0.8, font = 2)
}

## Specificity
plot.new() # to handle the missing 8. plot
plot.new()
text(0.5, 0.5, "No change", cex = 1.4, font = 2)

# Loop through each method to plot
for (i in 1:length(alg_labels)) {
  alg <- unname(alg_labels[i])
  temp <- radar_lognorm_samp_bias %>%
    filter(algorithm == alg & change == "No change")
  radar_lognorm_alg <- rbind(rep(1, 4), rep(0, 4), temp[, c(-1, -2, -3)])
  rownames(radar_lognorm_alg) <- c("Max", "Min", as.character(temp$nobs_cat))
  
  fmsb::radarchart(
    radar_lognorm_alg,
    axistype = 1,
    vlabels = est_labels_short,
    pcol = nobs_col,
    # pfcol = alpha(nobs_col, 0.2),
    plty = 1:5,
    plwd = 2,
    cglcol = "grey",
    cglty = 1,
    axislabcol = "black",
    caxislabels = c("0", "", "0.5", "", "1"),
    calcex = 0.9,
    vlcex = 0.9,
    title = "",
  )
  
  # Add method title
  mtext(unname(alg_labels[i]), side = 3, line = -0.3, cex = 0.8, font = 2)
}

# Add legend for sample sizes
plot.new()
legend("center",
       legend = unique(radar_lognorm_samp_bias$nobs_cat),
       col = nobs_col,
       lty = 1:5,
       lwd = 2,
       bty = "n",
       cex = 0.8,
       title = "Sample size",
       title.font = 2
)
dev.off()


#* Figure S2 --------------------------------------------------------------------
# Radar plots of scaled and inverted MSE for quantifying systematic changes in
# log-normally distributed data.
layout_matrix <- matrix(c(
  rep(1, 4), # row 1: Sensitivity title
  2:5, # row 2: radar plots 1-4 Sensitivity
  6:9, # row 3: radar plots 5-7 Sensitivity
  rep(10, 4), # row 4: Specificity title
  11:14, # row 5: radar plots 1-4 Specificity
  15:18 # row 6: radar plots 5-7 Specificity
), nrow = 6, byrow = TRUE)
row_heights <- c(0.2, 1, 1, 0.2, 1, 1)

radar_lognorm_samp <- res %>%
  filter(distribution == "lognorm") %>%
  group_by(change, algorithm, nobs_cat) %>%
  summarise(
    Bias_range = mean(range_dif/sd_y, na.rm = TRUE),
    MSE_range = mean(range_dif^2/sd_y^2, na.rm = TRUE),
    Bias_variance = mean(var_dif/sd_y, na.rm = TRUE),
    MSE_variance = mean(var_dif^2/sd_y^2, na.rm = TRUE),
    Bias_MADM = mean(madm_dif/sd_y, na.rm = TRUE),
    MSE_MADM = mean(madm_dif^2/sd_y^2, na.rm = TRUE),
    Bias_NCP = mean(n_cpts_dif, na.rm = TRUE),
    MSE_NCP = mean(n_cpts_dif^2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # abs() is needed for bias and does not do harm for MSE
  mutate(
    across(
      matches("^(Bias|MSE)"),
      ~ 1 - scales::rescale(
        abs(.x),
        from = c(0, max(abs(.x), na.rm = TRUE))
      )
    )
  )

radar_lognorm_samp_mse <- radar_lognorm_samp %>%
  select(change, algorithm, nobs_cat, starts_with("MSE_"))

pdf("results/figures/FigS2_radar_lognorm_samp_mse.pdf",
    width = 6.69, #unit: inches (equals 170mm)
    height = 8 #maximum height: 225mm (8.85 in) for figure and legend
)
# Set up plotting area
layout(layout_matrix, heights = row_heights)
par(mar = c(1, 1, 1, 1), oma = c(1, 1, 1, 1), xpd = TRUE)

## Sensitivity
plot.new()
text(0.5, 0.5, "Change", cex = 1.4, font = 2)

# Loop through each method to plot
for (i in 1:length(alg_labels)) {
  alg <- unname(alg_labels[i])
  temp <- radar_lognorm_samp_mse %>%
    filter(algorithm == alg & change == "Change")
  radar_lognorm_alg <- rbind(rep(1, 4), 
                          rep(0, 4), 
                          temp[, c(-1, -2, -3)])
  rownames(radar_lognorm_alg) <- c("Max", "Min", as.character(temp$nobs_cat))
  
  fmsb::radarchart(
    radar_lognorm_alg,
    axistype = 1,
    vlabels = est_labels_short,
    pcol = nobs_col,
    plty = 1:5,
    plwd = 2,
    cglcol = "grey",
    cglty = 1,
    axislabcol = "black",
    caxislabels = c("0", "", "0.5", "", "1"),
    calcex = 0.9,
    vlcex = 0.9,
    title = "",
  )
  
  # Add method title
  mtext(unname(alg_labels[i]), side = 3, line = -0.3, cex = 0.8, font = 2)
}

## Specificity
plot.new() # to handle the missing 8. plot
plot.new()
text(0.5, 0.5, "No change", cex = 1.4, font = 2)

# Loop through each method to plot
for (i in 1:length(alg_labels)) {
  alg <- unname(alg_labels[i])
  temp <- radar_lognorm_samp_mse %>%
    filter(algorithm == alg & change == "No change")
  radar_lognorm_alg <- rbind(rep(1, 4), rep(0, 4), temp[, c(-1, -2, -3)])
  rownames(radar_lognorm_alg) <- c("Max", "Min", as.character(temp$nobs_cat))
  
  fmsb::radarchart(
    radar_lognorm_alg,
    axistype = 1,
    vlabels = est_labels_short,
    pcol = nobs_col,
    # pfcol = alpha(nobs_col, 0.2),
    plty = 1:5,
    plwd = 2,
    cglcol = "grey",
    cglty = 1,
    axislabcol = "black",
    caxislabels = c("0", "", "0.5", "", "1"),
    calcex = 0.9,
    vlcex = 0.9,
    title = "",
  )
  
  # Add method title
  mtext(unname(alg_labels[i]), side = 3, line = -0.3, cex = 0.8, font = 2)
}

# Add legend for sample sizes
plot.new()
legend("center",
       legend = unique(radar_lognorm_samp_mse$nobs_cat),
       col = nobs_col,
       lty = 1:5,
       lwd = 2,
       bty = "n",
       cex = 0.8,
       title = "Sample size",
       title.font = 2
)
dev.off()


#* Figure S3 -------------------------------------------------------------------
# Heatmap of the bias for range by systematic change pattern across sample sizes
# for normally distributed data.
rangedat <- res %>%
  dplyr::select(distribution, pattern, algorithm, nobs, snr, dmaxf, range_dif, sd_y) %>%
  mutate(nobs_cat = cut(nobs,
    breaks = seq(min(nobs), max(nobs) + 10, by = 10),
    right = FALSE,
    labels = as.character(seq(30, max(nobs), by = 10))
  )) %>%
  mutate(nobs_cat = as.numeric(as.character(nobs_cat))) %>%
  group_by(nobs_cat, distribution, pattern, algorithm) %>%
  summarise(mean_range_dif = mean(range_dif/sd_y), .groups = "drop") # calculate group mean

ggplot(rangedat[rangedat$distribution == "norm", ], aes(y = nobs_cat, x = algorithm, fill = mean_range_dif)) +
  geom_tile(color = "white") +
  scale_fill_diverging_fixed0(
    data = rangedat[rangedat$distribution == "norm", ],
    var = mean_range_dif,
    name = "Bias"
  ) +
  facet_wrap(~pattern, scales = "free_x", ncol = 5) +
  scale_y_continuous(
    breaks = seq(100, max(rangedat$nobs_cat), by = 100),
    minor_breaks = seq(50, max(rangedat$nobs_cat), by = 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    expand = expansion(add = 0) # Removes padding before first and after last category
  ) +
  labs(y = "Sample size", x = "Algorithm") + # , title = "Normal distribution: Bias range") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(),
    axis.text.y = element_text(),
    axis.ticks.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )
ggsave("results/figures/FigS3_heatmap_bias_norm_range.pdf", width = 170, height = 200, units = "mm", bg = "white")


#* Figure S4 -------------------------------------------------------------------
# Heatmap of the bias for variance by systematic change pattern across sample sizes
# for normally distributed data.
vardat <- res %>%
  dplyr::select(distribution, pattern, algorithm, nobs, snr, dmaxf, var_dif, sd_y) %>%
  mutate(nobs_cat = cut(nobs,
    breaks = seq(min(nobs), max(nobs) + 10, by = 10),
    right = FALSE,
    labels = as.character(seq(30, max(nobs), by = 10))
  )) %>%
  mutate(nobs_cat = as.numeric(as.character(nobs_cat))) %>%
  group_by(nobs_cat, distribution, pattern, algorithm) %>%
  summarise(mean_var_dif = mean(var_dif/sd_y), .groups = "drop") # calculate group mean

ggplot(vardat[vardat$distribution == "norm", ], aes(y = nobs_cat, x = algorithm, fill = mean_var_dif)) +
  geom_tile(color = "white") +
  scale_fill_diverging_fixed0(
    data = vardat[vardat$distribution == "norm", ],
    var = mean_var_dif,
    name = "Bias"
  ) +
  facet_wrap(~pattern, scales = "free_x", ncol = 5) +
  scale_y_continuous(
    breaks = seq(100, max(vardat$nobs_cat), by = 100),
    minor_breaks = seq(50, max(vardat$nobs_cat), by = 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    expand = expansion(add = 0) # Removes padding before first and after last category
  ) +
  labs(y = "Sample size", x = "Algorithm") + # , title = "Normal distribution: Bias variance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(),
    axis.text.y = element_text(),
    axis.ticks.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )
ggsave("results/figures/FigS4_heatmap_bias_norm_var.pdf", width = 170, height = 200, units = "mm", bg = "white")


#* Figure S5 -------------------------------------------------------------------
# Heatmap of the bias for mean absolute deviation around the median by systematic
# change pattern across sample sizes for normally distributed data.
madmdat <- res %>%
  dplyr::select(distribution, pattern, algorithm, nobs, snr, dmaxf, madm_dif, sd_y) %>%
  mutate(nobs_cat = cut(nobs,
    breaks = seq(min(nobs), max(nobs) + 10, by = 10),
    right = FALSE,
    labels = as.character(seq(30, max(nobs), by = 10))
  )) %>%
  mutate(nobs_cat = as.numeric(as.character(nobs_cat))) %>%
  group_by(nobs_cat, distribution, pattern, algorithm) %>%
  summarise(mean_madm_dif = mean(madm_dif/sd_y), .groups = "drop") # calculate group mean

ggplot(madmdat[madmdat$distribution == "norm", ], aes(y = nobs_cat, x = algorithm, fill = mean_madm_dif)) +
  geom_tile(color = "white") +
  scale_fill_diverging_fixed0(
    data = madmdat[madmdat$distribution == "norm", ],
    var = mean_madm_dif,
    name = "Bias"
  ) +
  facet_wrap(~pattern, scales = "free_x", ncol = 5) +
  scale_y_continuous(
    breaks = seq(100, max(madmdat$nobs_cat), by = 100),
    minor_breaks = seq(50, max(madmdat$nobs_cat), by = 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    expand = expansion(add = 0) # Removes padding before first and after last category
  ) +
  labs(y = "Sample size", x = "Algorithm") + # , title = "Normal distribution: Bias MADM") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(),
    axis.text.y = element_text(),
    axis.ticks.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )
ggsave("results/figures/FigS5_heatmap_bias_norm_madm.pdf", width = 170, height = 200, units = "mm", bg = "white")


#* Figure S6 -------------------------------------------------------------------
# Heatmap of the bias for number of change points by systematic change pattern
# across sample sizes for normally distributed data.
ncptsdat <- res %>%
  dplyr::select(distribution, pattern, algorithm, nobs, snr, dmaxf, n_cpts_dif) %>%
  mutate(nobs_cat = cut(nobs,
    breaks = seq(min(nobs), max(nobs) + 10, by = 10),
    right = FALSE,
    labels = as.character(seq(30, max(nobs), by = 10))
  )) %>%
  mutate(nobs_cat = as.numeric(as.character(nobs_cat))) %>%
  group_by(nobs_cat, distribution, pattern, algorithm) %>%
  summarise(mean_n_cpts_dif = mean(n_cpts_dif), .groups = "drop") # calculate group mean

ggplot(ncptsdat[ncptsdat$distribution == "norm", ], aes(y = nobs_cat, x = algorithm, fill = mean_n_cpts_dif)) +
  geom_tile(color = "white") +
  scale_fill_diverging_fixed0(
    data = ncptsdat[ncptsdat$distribution == "norm", ],
    var = mean_n_cpts_dif,
    name = "Bias"
  ) +
  facet_wrap(~pattern, scales = "free_x", ncol = 5) +
  scale_y_continuous(
    breaks = seq(100, max(ncptsdat$nobs_cat), by = 100),
    minor_breaks = seq(50, max(ncptsdat$nobs_cat), by = 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    expand = expansion(add = 0) # Removes padding before first and after last category
  ) +
  labs(y = "Sample size", x = "Algorithm") + # , title = "Normal distribution: Bias NCP") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(),
    axis.text.y = element_text(),
    axis.ticks.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )
ggsave("results/figures/FigS6_heatmap_bias_norm_ncpts.pdf", width = 170, height = 200, units = "mm", bg = "white")


#* Figure S7 -------------------------------------------------------------------
# Heatmap of the mean squared error for range by systematic change pattern across
# sample sizes for normally distributed data.
rangedat <- res %>%
  dplyr::select(distribution, pattern, algorithm, nobs, snr, dmaxf, range_dif, sd_y) %>%
  mutate(nobs_cat = cut(nobs,
    breaks = seq(min(nobs), max(nobs) + 10, by = 10),
    right = FALSE,
    labels = as.character(seq(30, max(nobs), by = 10))
  )) %>%
  mutate(nobs_cat = as.numeric(as.character(nobs_cat))) %>%
  group_by(nobs_cat, distribution, pattern, algorithm) %>%
  summarise(mean_range_dif = mean(range_dif^2/sd_y^2), .groups = "drop") # calculate group mean

ggplot(rangedat[rangedat$distribution == "norm", ], aes(y = nobs_cat, x = algorithm, fill = mean_range_dif)) +
  geom_tile(color = "white") +
  scale_fill_positive_fixed0(
    data = rangedat[rangedat$distribution == "norm", ],
    var = mean_range_dif,
    name = "MSE"
  ) +
  facet_wrap(~pattern, scales = "free_x", ncol = 5) +
  scale_y_continuous(
    breaks = seq(100, max(rangedat$nobs_cat), by = 100),
    minor_breaks = seq(50, max(rangedat$nobs_cat), by = 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    expand = expansion(add = 0) # Removes padding before first and after last category
  ) +
  labs(y = "Sample size", x = "Algorithm") + # , title = "Normal distribution: MSE range") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(),
    axis.text.y = element_text(),
    axis.ticks.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )
ggsave("results/figures/FigS7_heatmap_mse_norm_range.pdf", width = 170, height = 200, units = "mm", bg = "white")


#* Figure S8 -------------------------------------------------------------------
# Heatmap of the mean squared error for variance by systematic change pattern
# across sample sizes for normally distributed data.
vardat <- res %>%
  dplyr::select(distribution, pattern, algorithm, nobs, snr, dmaxf, var_dif, sd_y) %>%
  mutate(nobs_cat = cut(nobs,
    breaks = seq(min(nobs), max(nobs) + 10, by = 10),
    right = FALSE,
    labels = as.character(seq(30, max(nobs), by = 10))
  )) %>%
  mutate(nobs_cat = as.numeric(as.character(nobs_cat))) %>%
  group_by(nobs_cat, distribution, pattern, algorithm) %>%
  summarise(mean_var_dif = mean(var_dif^2/sd_y^2), .groups = "drop") # calculate group mean

ggplot(vardat[vardat$distribution == "norm", ], aes(y = nobs_cat, x = algorithm, fill = mean_var_dif)) +
  geom_tile(color = "white") +
  scale_fill_positive_fixed0(
    data = vardat[vardat$distribution == "norm", ],
    var = mean_var_dif,
    name = "MSE"
  ) +
  facet_wrap(~pattern, scales = "free_x", ncol = 5) +
  scale_y_continuous(
    breaks = seq(100, max(vardat$nobs_cat), by = 100),
    minor_breaks = seq(50, max(vardat$nobs_cat), by = 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    expand = expansion(add = 0) # Removes padding before first and after last category
  ) +
  labs(y = "Sample size", x = "Algorithm") + # , title = "Normal distribution: MSE variance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(),
    axis.text.y = element_text(),
    axis.ticks.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )
ggsave("results/figures/FigS8_heatmap_mse_norm_var.pdf", width = 170, height = 200, units = "mm", bg = "white")


#* Figure S9 -------------------------------------------------------------------
# Heatmap of the mean squared error for mean absolute deviation around the median
# by systematic change pattern across sample sizes for normally distributed data.
madmdat <- res %>%
  dplyr::select(distribution, pattern, algorithm, nobs, snr, dmaxf, madm_dif, sd_y) %>%
  mutate(nobs_cat = cut(nobs,
    breaks = seq(min(nobs), max(nobs) + 10, by = 10),
    right = FALSE,
    labels = as.character(seq(30, max(nobs), by = 10))
  )) %>%
  mutate(nobs_cat = as.numeric(as.character(nobs_cat))) %>%
  group_by(nobs_cat, distribution, pattern, algorithm) %>%
  summarise(mean_madm_dif = mean(madm_dif^2/sd_y^2), .groups = "drop") # calculate group mean

ggplot(madmdat[madmdat$distribution == "norm", ], aes(y = nobs_cat, x = algorithm, fill = mean_madm_dif)) +
  geom_tile(color = "white") +
  scale_fill_positive_fixed0(
    data = madmdat[madmdat$distribution == "norm", ],
    var = mean_madm_dif,
    name = "MSE"
  ) +
  facet_wrap(~pattern, scales = "free_x", ncol = 5) +
  scale_y_continuous(
    breaks = seq(100, max(madmdat$nobs_cat), by = 100),
    minor_breaks = seq(50, max(madmdat$nobs_cat), by = 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    expand = expansion(add = 0) # Removes padding before first and after last category
  ) +
  labs(y = "Sample size", x = "Algorithm") + # , title = "Normal distribution: MSE MADM") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(),
    axis.text.y = element_text(),
    axis.ticks.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )
ggsave("results/figures/FigS9_heatmap_mse_norm_madm.pdf", width = 170, height = 200, units = "mm", bg = "white")


#* Figure S10 -------------------------------------------------------------------
# Heatmap of the mean squared error for number of change points by systematic
# change pattern across sample sizes for normally distributed data.
ncptsdat <- res %>%
  dplyr::select(distribution, pattern, algorithm, nobs, snr, dmaxf, n_cpts_dif) %>%
  mutate(nobs_cat = cut(nobs,
    breaks = seq(min(nobs), max(nobs) + 10, by = 10),
    right = FALSE,
    labels = as.character(seq(30, max(nobs), by = 10))
  )) %>%
  mutate(nobs_cat = as.numeric(as.character(nobs_cat))) %>%
  group_by(nobs_cat, distribution, pattern, algorithm) %>%
  summarise(mean_n_cpts_dif = mean(n_cpts_dif^2), .groups = "drop") # calculate group mean

ggplot(ncptsdat[ncptsdat$distribution == "norm", ], aes(y = nobs_cat, x = algorithm, fill = mean_n_cpts_dif)) +
  geom_tile(color = "white") +
  scale_fill_positive_fixed0(
    data = ncptsdat[ncptsdat$distribution == "norm", ],
    var = mean_n_cpts_dif,
    name = "MSE"
  ) +
  facet_wrap(~pattern, scales = "free_x", ncol = 5) +
  scale_y_continuous(
    breaks = seq(100, max(ncptsdat$nobs_cat), by = 100),
    minor_breaks = seq(50, max(ncptsdat$nobs_cat), by = 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    expand = expansion(add = 0) # Removes padding before first and after last category
  ) +
  labs(y = "Sample size", x = "Algorithm") + # , title = "Normal distribution: MSE NCP") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(),
    axis.text.y = element_text(),
    axis.ticks.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )
ggsave("results/figures/FigS10_heatmap_mse_norm_ncpts.pdf", width = 170, height = 200, units = "mm", bg = "white")


#* Figure S11 -------------------------------------------------------------------
# Heatmap of the bias for range by systematic change pattern across sample sizes
# for log-normally distributed data.
rangedat <- res %>%
  dplyr::select(distribution, pattern, algorithm, nobs, snr, dmaxf, range_dif,sd_y) %>%
  mutate(nobs_cat = cut(nobs,
    breaks = seq(min(nobs), max(nobs) + 10, by = 10),
    right = FALSE,
    labels = as.character(seq(30, max(nobs), by = 10))
  )) %>%
  mutate(nobs_cat = as.numeric(as.character(nobs_cat))) %>%
  group_by(nobs_cat, distribution, pattern, algorithm) %>%
  summarise(mean_range_dif = mean(range_dif/sd_y), .groups = "drop") # calculate group mean

ggplot(rangedat[rangedat$distribution == "lognorm", ], aes(y = nobs_cat, x = algorithm, fill = mean_range_dif)) +
  geom_tile(color = "white") +
  scale_fill_diverging_fixed0(
    data = rangedat[rangedat$distribution == "lognorm", ],
    var = mean_range_dif,
    name = "Bias"
  ) +
  facet_wrap(~pattern, scales = "free_x", ncol = 5) +
  scale_y_continuous(
    breaks = seq(100, max(rangedat$nobs_cat), by = 100),
    minor_breaks = seq(50, max(rangedat$nobs_cat), by = 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    expand = expansion(add = 0) # Removes padding before first and after last category
  ) +
  labs(y = "Sample size", x = "Algorithm") + # , title = "Lognormal distribution: Bias range") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(),
    axis.text.y = element_text(),
    axis.ticks.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )
ggsave("results/figures/FigS10_heatmap_bias_lognorm_range.pdf", width = 170, height = 200, units = "mm", bg = "white")


#* Figure S12 ------------------------------------------------------------------
# Heatmap of the bias for variance by systematic change pattern across sample sizes
# for log-normally distributed data.
vardat <- res %>%
  dplyr::select(distribution, pattern, algorithm, nobs, snr, dmaxf, var_dif, sd_y) %>%
  mutate(nobs_cat = cut(nobs,
    breaks = seq(min(nobs), max(nobs) + 10, by = 10),
    right = FALSE,
    labels = as.character(seq(30, max(nobs), by = 10))
  )) %>%
  mutate(nobs_cat = as.numeric(as.character(nobs_cat))) %>%
  group_by(nobs_cat, distribution, pattern, algorithm) %>%
  summarise(mean_var_dif = mean(var_dif/sd_y), .groups = "drop") # calculate group mean

ggplot(vardat[vardat$distribution == "lognorm", ], aes(y = nobs_cat, x = algorithm, fill = mean_var_dif)) +
  geom_tile(color = "white") +
  scale_fill_diverging_fixed0(
    data = vardat[vardat$distribution == "lognorm", ],
    var = mean_var_dif,
    name = "Bias"
  ) +
  facet_wrap(~pattern, scales = "free_x", ncol = 5) +
  scale_y_continuous(
    breaks = seq(100, max(vardat$nobs_cat), by = 100),
    minor_breaks = seq(50, max(vardat$nobs_cat), by = 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    expand = expansion(add = 0) # Removes padding before first and after last category
  ) +
  labs(y = "Sample size", x = "Algorithm") + # , title = "Lognormal distribution: Bias variance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(),
    axis.text.y = element_text(),
    axis.ticks.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )
ggsave("results/figures/FigS12_heatmap_bias_lognorm_var.pdf", width = 170, height = 200, units = "mm", bg = "white")


#* Figure S13 ------------------------------------------------------------------
# Heatmap of the bias for mean absolute deviation around the median by systematic
# change pattern across sample sizes for log-normally distributed data.
madmdat <- res %>%
  dplyr::select(distribution, pattern, algorithm, nobs, snr, dmaxf, madm_dif, sd_y) %>%
  mutate(nobs_cat = cut(nobs,
    breaks = seq(min(nobs), max(nobs) + 10, by = 10),
    right = FALSE,
    labels = as.character(seq(30, max(nobs), by = 10))
  )) %>%
  mutate(nobs_cat = as.numeric(as.character(nobs_cat))) %>%
  group_by(nobs_cat, distribution, pattern, algorithm) %>%
  summarise(mean_madm_dif = mean(madm_dif/sd_y), .groups = "drop") # calculate group mean

ggplot(madmdat[madmdat$distribution == "lognorm", ], aes(y = nobs_cat, x = algorithm, fill = mean_madm_dif)) +
  geom_tile(color = "white") +
  scale_fill_diverging_fixed0(
    data = madmdat[madmdat$distribution == "lognorm", ],
    var = mean_madm_dif,
    name = "Bias"
  ) +
  facet_wrap(~pattern, scales = "free_x", ncol = 5) +
  scale_y_continuous(
    breaks = seq(100, max(madmdat$nobs_cat), by = 100),
    minor_breaks = seq(50, max(madmdat$nobs_cat), by = 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    expand = expansion(add = 0) # Removes padding before first and after last category
  ) +
  labs(y = "Sample size", x = "Algorithm") + # , title = "Lognormal distribution: Bias MADM") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(),
    axis.text.y = element_text(),
    axis.ticks.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )
ggsave("results/figures/FigS13_heatmap_bias_lognorm_madm.pdf", width = 170, height = 200, units = "mm", bg = "white")


#* Figure S14 ------------------------------------------------------------------
# Heatmap of the bias for number of change points by systematic change pattern
# across sample sizes for log-normally distributed data.
ncptsdat <- res %>%
  dplyr::select(distribution, pattern, algorithm, nobs, snr, dmaxf, n_cpts_dif) %>%
  mutate(nobs_cat = cut(nobs,
    breaks = seq(min(nobs), max(nobs) + 10, by = 10),
    right = FALSE,
    labels = as.character(seq(30, max(nobs), by = 10))
  )) %>%
  mutate(nobs_cat = as.numeric(as.character(nobs_cat))) %>%
  group_by(nobs_cat, distribution, pattern, algorithm) %>%
  summarise(mean_n_cpts_dif = mean(n_cpts_dif), .groups = "drop") # calculate group mean

ggplot(ncptsdat[ncptsdat$distribution == "lognorm", ], aes(y = nobs_cat, x = algorithm, fill = mean_n_cpts_dif)) +
  geom_tile(color = "white") +
  scale_fill_diverging_fixed0(
    data = ncptsdat[ncptsdat$distribution == "lognorm", ],
    var = mean_n_cpts_dif,
    name = "Bias"
  ) +
  facet_wrap(~pattern, scales = "free_x", ncol = 5) +
  scale_y_continuous(
    breaks = seq(100, max(ncptsdat$nobs_cat), by = 100),
    minor_breaks = seq(50, max(ncptsdat$nobs_cat), by = 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    expand = expansion(add = 0) # Removes padding before first and after last category
  ) +
  labs(y = "Sample size", x = "Algorithm") + # , title = "Lognormal distribution: Bias NCP") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(),
    axis.text.y = element_text(),
    axis.ticks.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )
ggsave("results/figures/FigS14_heatmap_bias_lognorm_ncpts.pdf", width = 170, height = 200, units = "mm", bg = "white")


#* Figure S15 ------------------------------------------------------------------
# Heatmap of the mean squared error for range by systematic change pattern across
# sample sizes for log-normally distributed data.
rangedat <- res %>%
  dplyr::select(distribution, pattern, algorithm, nobs, snr, dmaxf, range_dif, sd_y) %>%
  mutate(nobs_cat = cut(nobs,
    breaks = seq(min(nobs), max(nobs) + 10, by = 10),
    right = FALSE,
    labels = as.character(seq(30, max(nobs), by = 10))
  )) %>%
  mutate(nobs_cat = as.numeric(as.character(nobs_cat))) %>%
  group_by(nobs_cat, distribution, pattern, algorithm) %>%
  summarise(mean_range_dif = mean(range_dif^2/sd_y^2), .groups = "drop") # calculate group mean

ggplot(rangedat[rangedat$distribution == "lognorm", ], aes(y = nobs_cat, x = algorithm, fill = mean_range_dif)) +
  geom_tile(color = "white") +
  scale_fill_positive_fixed0(
    data = rangedat[rangedat$distribution == "lognorm", ],
    var = mean_range_dif,
    name = "MSE"
  ) +
  facet_wrap(~pattern, scales = "free_x", ncol = 5) +
  scale_y_continuous(
    breaks = seq(100, max(rangedat$nobs_cat), by = 100),
    minor_breaks = seq(50, max(rangedat$nobs_cat), by = 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    expand = expansion(add = 0) # Removes padding before first and after last category
  ) +
  labs(y = "Sample size", x = "Algorithm") + # , title = "Lognormal distribution: MSE range") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(),
    axis.text.y = element_text(),
    axis.ticks.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )
ggsave("results/figures/FigS15_heatmap_mse_lognorm_range.pdf", width = 170, height = 200, units = "mm", bg = "white")


#* Figure S16 ------------------------------------------------------------------
# Heatmap of the mean squared error for variance by systematic change pattern
# across sample sizes for log-normally distributed data.
vardat <- res %>%
  dplyr::select(distribution, pattern, algorithm, nobs, snr, dmaxf, var_dif, sd_y) %>%
  mutate(nobs_cat = cut(nobs,
    breaks = seq(min(nobs), max(nobs) + 10, by = 10),
    right = FALSE,
    labels = as.character(seq(30, max(nobs), by = 10))
  )) %>%
  mutate(nobs_cat = as.numeric(as.character(nobs_cat))) %>%
  group_by(nobs_cat, distribution, pattern, algorithm) %>%
  summarise(mean_var_dif = mean(var_dif^2/sd_y^2), .groups = "drop") # calculate group mean

ggplot(vardat[vardat$distribution == "lognorm", ], aes(y = nobs_cat, x = algorithm, fill = mean_var_dif)) +
  geom_tile(color = "white") +
  scale_fill_positive_fixed0(
    data = vardat[vardat$distribution == "lognorm", ],
    var = mean_var_dif,
    name = "MSE"
  ) +
  facet_wrap(~pattern, scales = "free_x", ncol = 5) +
  scale_y_continuous(
    breaks = seq(100, max(vardat$nobs_cat), by = 100),
    minor_breaks = seq(50, max(vardat$nobs_cat), by = 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    expand = expansion(add = 0) # Removes padding before first and after last category
  ) +
  labs(y = "Sample size", x = "Algorithm") + # , title = "Lognormal distribution: MSE variance") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(),
    axis.text.y = element_text(),
    axis.ticks.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )
ggsave("results/figures/FigS16_heatmap_mse_lognorm_var.pdf", width = 170, height = 200, units = "mm", bg = "white")


#* Figure S17 ------------------------------------------------------------------
# Heatmap of the mean squared error for mean absolute deviation around the median
# by systematic change pattern across sample sizes for log-normally distributed data.
madmdat <- res %>%
  dplyr::select(distribution, pattern, algorithm, nobs, snr, dmaxf, madm_dif, sd_y) %>%
  mutate(nobs_cat = cut(nobs,
    breaks = seq(min(nobs), max(nobs) + 10, by = 10),
    right = FALSE,
    labels = as.character(seq(30, max(nobs), by = 10))
  )) %>%
  mutate(nobs_cat = as.numeric(as.character(nobs_cat))) %>%
  group_by(nobs_cat, distribution, pattern, algorithm) %>%
  summarise(mean_madm_dif = mean(madm_dif^2/sd_y^2), .groups = "drop") # calculate group mean

ggplot(madmdat[madmdat$distribution == "lognorm", ], aes(y = nobs_cat, x = algorithm, fill = mean_madm_dif)) +
  geom_tile(color = "white") +
  scale_fill_positive_fixed0(
    data = madmdat[madmdat$distribution == "lognorm", ],
    var = mean_madm_dif,
    name = "MSE"
  ) +
  facet_wrap(~pattern, scales = "free_x", ncol = 5) +
  scale_y_continuous(
    breaks = seq(100, max(madmdat$nobs_cat), by = 100),
    minor_breaks = seq(50, max(madmdat$nobs_cat), by = 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    expand = expansion(add = 0) # Removes padding before first and after last category
  ) +
  labs(y = "Sample size", x = "Algorithm") + # , title = "Lognormal distribution: MSE MADM") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(),
    axis.text.y = element_text(),
    axis.ticks.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )
ggsave("results/figures/FigS17_heatmap_mse_lognorm_madm.pdf", width = 170, height = 200, units = "mm", bg = "white")


#* Figure S18 ------------------------------------------------------------------
# Heatmap of the mean squared error for number of change points by systematic
# change pattern across sample sizes for log-normally distributed data.
ncptsdat <- res %>%
  dplyr::select(distribution, pattern, algorithm, nobs, snr, dmaxf, n_cpts_dif) %>%
  mutate(nobs_cat = cut(nobs,
    breaks = seq(min(nobs), max(nobs) + 10, by = 10),
    right = FALSE,
    labels = as.character(seq(30, max(nobs), by = 10))
  )) %>%
  mutate(nobs_cat = as.numeric(as.character(nobs_cat))) %>%
  group_by(nobs_cat, distribution, pattern, algorithm) %>%
  summarise(mean_n_cpts_dif = mean(n_cpts_dif^2), .groups = "drop") # calculate group mean

ggplot(ncptsdat[ncptsdat$distribution == "lognorm", ], aes(y = nobs_cat, x = algorithm, fill = mean_n_cpts_dif)) +
  geom_tile(color = "white") +
  scale_fill_positive_fixed0(
    data = ncptsdat[ncptsdat$distribution == "lognorm", ],
    var = mean_n_cpts_dif,
    name = "MSE"
  ) +
  facet_wrap(~pattern, scales = "free_x", ncol = 5) +
  scale_y_continuous(
    breaks = seq(100, max(ncptsdat$nobs_cat), by = 100),
    minor_breaks = seq(50, max(ncptsdat$nobs_cat), by = 50),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    expand = expansion(add = 0) # Removes padding before first and after last category
  ) +
  labs(y = "Sample size", x = "Algorithm") + # , title = "Lognormal distribution: MSE NCP") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.ticks.x = element_blank(),
    axis.title.x = element_text(),
    axis.text.y = element_text(),
    axis.ticks.y = element_line(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.spacing = unit(1, "lines"),
    strip.text = element_text(face = "bold")
  )
ggsave("results/figures/FigS16_heatmap_mse_lognorm_ncpts.pdf", width = 170, height = 200, units = "mm", bg = "white")