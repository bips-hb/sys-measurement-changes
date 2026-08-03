## Additional file 1
# excel file with detailed bias and mean squared error

# load packages & functions
library(dplyr)
library(openxlsx2)


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

# change pattern
pattern_labels <- c(
  "No change",
  "Jump",
  "Linear",
  "Quadratic",
  "Recurrent"
)

## preprocess results
res <- res %>%
  mutate(algorithm = factor(
    algorithm,
    levels = names(alg_labels),
    labels = unname(alg_labels)
  )) %>%
  mutate(pattern = factor(
    pattern,
    levels = c("nochange", "jump", "linear", "quadratic", "recurrent"),
    labels = pattern_labels
  ))


# Add. file 1 ------------------------------------------------------------------

res_norm_bias_mse <- res %>%
  filter(distribution == "norm") %>%
  group_by(nobs_cat, pattern, dmaxf, snr, algorithm) %>%
  summarise(
    N = n(),
    Bias_range = mean(range_dif),
    Bias_range_SE = sd(range_dif)/sqrt(n()),
    MSE_range = mean(range_dif^2),
    MSE_range_SE = sd(range_dif^2)/sqrt(n()),
    Bias_variance = mean(var_dif),
    Bias_variance_SE = sd(var_dif)/sqrt(n()),
    MSE_variance = mean(var_dif^2),
    MSE_variance_SE = sd(var_dif^2)/sqrt(n()),
    Bias_MADM = mean(madm_dif),
    Bias_MADM_SE = sd(madm_dif)/sqrt(n()),
    MSE_MADM = mean(madm_dif^2),
    MSE_MADM_SE = sd(madm_dif^2)/sqrt(n()),
    Bias_NCP = mean(n_cpts_dif),
    Bias_NCP_SE = sd(n_cpts_dif)/sqrt(n()),
    MSE_NCP = mean(n_cpts_dif^2),
    MSE_NCP_SE = sd(n_cpts_dif^2)/sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(across(
    matches("^(Bias|MSE)") & !ends_with("_SE"),
    ~ rank(abs(.x), ties.method = "average"),
    # abs() is needed for bias and does not do harm for MSE
    .names = "{.col}_rank"
  )) %>%
  mutate(Bias_borda = rowSums(across(matches("^Bias.*_rank$"))),
         Bias_borda_scale = (Bias_borda - 4*1) / (4*7 - 4*1), 
         MSE_borda = rowSums(across(matches("^MSE.*_rank$"))),
         MSE_borda_scale = (MSE_borda - 4*1) / (4*7 - 4*1) ) %>% 
  ungroup() %>%
  select(-c(ends_with("_rank"), "Bias_borda", "MSE_borda")) %>%
  rename(
    Sample_size = nobs_cat,
    Change_pattern = pattern,
    Magnitude_change = dmaxf,
    SNR = snr,
    Method = algorithm,
    Bias_Borda = Bias_borda_scale,
    MSE_Borda = MSE_borda_scale
  )

res_lognorm_bias_mse <- res %>%
  filter(distribution == "lognorm") %>%
  group_by(nobs_cat, pattern, dmaxf, snr, algorithm) %>%
  summarise(
    N = n(),
    Bias_range = mean(range_dif),
    Bias_range_SE = sd(range_dif)/sqrt(n()),
    MSE_range = mean(range_dif^2),
    MSE_range_SE = sd(range_dif^2)/sqrt(n()),
    Bias_variance = mean(var_dif),
    Bias_variance_SE = sd(var_dif)/sqrt(n()),
    MSE_variance = mean(var_dif^2),
    MSE_variance_SE = sd(var_dif^2)/sqrt(n()),
    Bias_MADM = mean(madm_dif),
    Bias_MADM_SE = sd(madm_dif)/sqrt(n()),
    MSE_MADM = mean(madm_dif^2),
    MSE_MADM_SE = sd(madm_dif^2)/sqrt(n()),
    Bias_NCP = mean(n_cpts_dif),
    Bias_NCP_SE = sd(n_cpts_dif)/sqrt(n()),
    MSE_NCP = mean(n_cpts_dif^2),
    MSE_NCP_SE = sd(n_cpts_dif^2)/sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(across(
    matches("^(Bias|MSE)") & !ends_with("_SE"),
    ~ rank(abs(.x), ties.method = "average"),  
    # abs() is needed for bias and does not do harm for MSE
    .names = "{.col}_rank"
  )) %>%
  mutate(Bias_borda = rowSums(across(matches("^Bias.*_rank$"))),
         Bias_borda_scale = (Bias_borda - 4*1) / (4*7 - 4*1), 
         MSE_borda = rowSums(across(matches("^MSE.*_rank$"))),
         MSE_borda_scale = (MSE_borda - 4*1) / (4*7 - 4*1) ) %>% 
  ungroup() %>%
  select(-c(ends_with("_rank"), "Bias_borda", "MSE_borda")) %>%
  rename(
    Sample_size = nobs_cat,
    Change_pattern = pattern,
    Magnitude_change = dmaxf,
    SNR = snr,
    Method = algorithm,
    Bias_Borda = Bias_borda_scale,
    MSE_Borda = MSE_borda_scale
  )

## save excel file
wb <- openxlsx2::wb_workbook()
border_col <- openxlsx2::wb_color(theme = 1)
border_sty <- "thin"

## sheet 1: Normal distribution
dims_header <- openxlsx2::wb_dims(
  rows = 1,
  cols = seq_len(ncol(res_norm_bias_mse))
)
dims_num <- openxlsx2::wb_dims(
  rows = 2:(nrow(res_norm_bias_mse) + 1),
  cols = seq(7, ncol(res_norm_bias_mse), by=1)
)
wb <- openxlsx2::wb_add_worksheet(wb, "Normal Distribution") %>%
  openxlsx2::wb_add_data(
    sheet = "Normal Distribution",
    x = res_norm_bias_mse
  ) %>%
  openxlsx2::wb_add_numfmt(
    dims = dims_num,
    numfmt = "0.000"
  ) %>%
  openxlsx2::wb_add_fill(
    dims = dims_header,
    color = openxlsx2::wb_color(hex = "#001260")
  ) %>%
  openxlsx2::wb_add_font(
    dims = dims_header,
    name = "Arial",
    bold = TRUE,
    color = openxlsx2::wb_color("white")
  )

# Auto-adjust column widths
wb <- openxlsx2::wb_set_col_widths(
  wb,
  sheet = "Normal Distribution",
  cols = 1:ncol(res_norm_bias_mse), widths = "auto"
)

## sheet 2: Lognormal distribution
dims_header <- openxlsx2::wb_dims(
  rows = 1,
  cols = seq_len(ncol(res_lognorm_bias_mse))
)
wb <- openxlsx2::wb_add_worksheet(wb, "Log-normal Distribution") %>%
  openxlsx2::wb_add_data(
    sheet = "Log-normal Distribution",
    x = res_lognorm_bias_mse
  ) %>%
  openxlsx2::wb_add_fill(
    dims = dims_header,
    color = openxlsx2::wb_color(hex = "#001260")
  ) %>%
  openxlsx2::wb_add_font(
    dims = dims_header,
    name = "Arial",
    bold = TRUE,
    color = openxlsx2::wb_color("white")
  )

# Auto-adjust column widths
wb <- openxlsx2::wb_set_col_widths(
  wb,
  sheet = "Log-normal Distribution",
  cols = 1:ncol(res_lognorm_bias_mse), widths = "auto"
)

openxlsx2::wb_save(wb, "results/additional_file_1.xlsx")

