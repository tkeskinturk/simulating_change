
# - promises and pitfalls of panel change ------------------------------------ #
# - ks distance example ------------------------------------------------------ #

### note: this script
###       (a) generates an exemplary panel data from a DGP A,
###       (b) simulates an exemplary panel data from a DGP B,
###       (c) showcases the KS distance approach using their slope distribution.

rm(list = ls()) # clean-up

# install.packages("devtools")
# devtools::install_github("tkeskinturk/gridsearch")

pacman::p_load(
  purrr,
  gridsearch,
  tidyverse,
  hrbrthemes,
  GGally,
  patchwork)
theme_set(theme_ipsum_rc())

# --- helpers
my_oka <-
  c("#eebd64",
    "#8b3b36", 
    "#a2b4b7")
`%nin%` = Negate(`%in%`)

set.seed(112358)

# --- glm function
glm_fit <- function(data) {
    fit <- suppressWarnings(
      stats::glm(y ~ scales::rescale(t), 
                 data = data, family = stats::binomial()))
    prob1 <- stats::plogis(fit$coefficients[1])
    prob2 <- stats::plogis(fit$coefficients[1] + fit$coefficients[2])
    estimate <- round(prob2 - prob1, 3)
    return(estimate)
  }

# PART 1: DGPs --------------------------------------------------------------- #

### DGP A

# --- read in
data_A <- buildDGP(
  n = 500,
  t = 12,
  rate = 0.5,
  balance_dir = 1,
  balance_res = .50,
  strength =  1,
  reliable = .8,
  export = TRUE
)$data

# --- organize
data_A <- data_A |>
  as_tibble() |>
  mutate(pid = 1:500) |>
  pivot_longer(cols = -pid,
               names_to = "cols",
               values_to = "y") |>
  mutate(cols = factor(cols, levels = paste0("V", 1:12))) |> 
  group_by(cols) |> mutate(t = cur_group_id()) |> ungroup() |>
  nest(.by = "pid") |>
  mutate(slopes_A = purrr::map_dbl(.x = data, .f = ~ glm_fit(.x))) |>
  select(pid, slopes_A)

### DGP B

# --- read in
data_B <- buildDGP(
  n = 500,
  t = 12,
  rate = 0.5,
  balance_dir = .50,
  balance_res = .50,
  strength = 1,
  reliable = .8,
  export = TRUE
)$data

# --- organize
data_B <- data_B |>
  as_tibble() |>
  mutate(pid = 1:500) |>
  pivot_longer(cols = -pid,
               names_to = "cols",
               values_to = "y") |>
  mutate(cols = factor(cols, levels = paste0("V", 1:12))) |> 
  group_by(cols) |> mutate(t = cur_group_id()) |> ungroup() |>
  nest(.by = "pid") |>
  mutate(slopes_B = purrr::map_dbl(.x = data, .f = ~ glm_fit(.x))) |>
  select(pid, slopes_B)

# PART 2: DISTRIBUTIONS ------------------------------------------------------ #

# --- plot 1
example_p1 <- data_A |>
  ggplot(aes(x = slopes_A, y = after_stat(ndensity))) +
  geom_histogram(
    bins = 15,
    color = "black",
    fill = my_oka[1],
    alpha = 0.8
  ) +
  theme_ipsum_rc(grid = "XxY", axis_title_size = 12) +
  labs(subtitle = "Slopes from DGP 1", 
       x = "", 
       y = "Relative Frequency") +
  scale_x_continuous(
    limits = c(-1.15, 1.15),
    oob = scales::oob_keep,
    breaks = scales::pretty_breaks()
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks())

# --- plot 2
example_p2 <- data_B |>
  ggplot(aes(x = slopes_B, y = after_stat(ndensity))) +
  geom_histogram(
    bins = 15,
    color = "black",
    fill = my_oka[3],
    alpha = 0.8
  ) +
  theme_ipsum_rc(grid = "XxY", axis_title_size = 12) +
  labs(subtitle = "Slopes from DGP 2", 
       x = "Individual-Level Slope Coefficients", 
       y = "") +
  scale_x_continuous(
    limits = c(-1.15, 1.15),
    oob = scales::oob_keep,
    breaks = scales::pretty_breaks()
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks())

# --- plot 3
example_p3 <-
  tibble(dgp1 = data_A$slopes_A, dgp2 = data_B$slopes_B) |> 
  ggplot() +
  stat_ecdf(aes(x = dgp1), geom = "step", col = my_oka[1]) +
  stat_ecdf(aes(x = dgp2), geom = "step", col = my_oka[3]) +
  theme_ipsum_rc(grid = "XxY", axis_title_size = 12) +
  labs(subtitle = "Empirical CDF",
       x = "Individual-Level Slope Coefficients", 
       y = "ECDF") +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks())

ks.test(data_A$slopes_A, data_B$slopes_B, exact = TRUE) # 0.164

# --- draw trajectories
png(
  "figures/ks_distance.png",
  w = 12,
  h =  7,
  units = "in",
  res = 500
)
(example_p1 / example_p2) | example_p3
dev.off()

# ---------------------------------------------------------------------------- #
