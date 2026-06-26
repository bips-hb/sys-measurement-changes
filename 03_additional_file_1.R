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
  alg_gam = "GAM",
  alg_lowess = "LOWESS",
  alg_pelt = "PELT",
  alg_piecewise_reg = "PR"
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
    Bias_range = mean(range_dif, na.rm = TRUE),
    MSE_range = mean(range_dif^2, na.rm = TRUE),
    Bias_variance = mean(var_dif, na.rm = TRUE),
    MSE_variance = mean(var_dif^2, na.rm = TRUE),
    Bias_MADM = mean(madm_dif, na.rm = TRUE),
    MSE_MADM = mean(madm_dif^2, na.rm = TRUE),
    Bias_NCP = mean(n_cpts_dif, na.rm = TRUE),
    MSE_NCP = mean(n_cpts_dif^2, na.rm = TRUE)
  ) %>%
  mutate(across(
    matches("^(Bias|MSE)"),
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
    Method = algorithm
  )

res_lognorm_bias_mse <- res %>%
  filter(distribution == "lognorm") %>%
  group_by(nobs_cat, pattern, dmaxf, snr, algorithm) %>%
  summarise(
    Bias_range = mean(range_dif, na.rm = TRUE),
    MSE_range = mean(range_dif^2, na.rm = TRUE),
    Bias_variance = mean(var_dif, na.rm = TRUE),
    MSE_variance = mean(var_dif^2, na.rm = TRUE),
    Bias_MADM = mean(madm_dif, na.rm = TRUE),
    MSE_MADM = mean(madm_dif^2, na.rm = TRUE),
    Bias_NCP = mean(n_cpts_dif, na.rm = TRUE),
    MSE_NCP = mean(n_cpts_dif^2, na.rm = TRUE)
  ) %>%
  mutate(across(
    matches("^(Bias|MSE)"),
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
    Method = algorithm
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
wb <- openxlsx2::wb_add_worksheet(wb, "Normal Distribution") %>%
  openxlsx2::wb_add_data(
    sheet = "Normal Distribution",
    x = res_norm_bias_mse
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

