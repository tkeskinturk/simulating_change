
# - promises and pitfalls of panel change ------------------------------------ #
# - case 4 # immigration ----------------------------------------------------- #

### note: this script
###       (a) takes `imm` as an examplary panel data item,
###       (b) simulates potential DGP distributions using fixed-effects models,
###       (c) shows the relative fit of the DGPs to observed data.

rm(list = ls()) # clean-up

# install.packages("devtools")
# devtools::install_github("tkeskinturk/gridsearch")

pacman::p_load(
  purrr,
  furrr,
  srvyr,
  hrbrthemes,
  gridsearch,
  tidyverse,
  patchwork)
theme_set(theme_ipsum_rc())

# --- helpers
my_oka <-
  c("#eebd64",
    "#8b3b36", 
    "#a2b4b7")
`%nin%` = Negate(`%in%`)

# PART 1: 26-WAVE DATA ------------------------------------------------------- #

# --- read in
data <- read_csv("./data/yllanon.csv")

# --- cleaning up
data <- data |>
  ## organize
  mutate(pid = cur_group_id(), .by = "id") |>
  ## chop off post-COVID periods
  filter(wave < 18) |>
  ## simplify
  select(pid, t = wave, y = imm) |>
  ## clean-up the outcome
  filter(y %nin% c(8, 9)) |> drop_na() |>
  ## fully balanced
  mutate(n = n(), .by = "pid") |> filter(n == 17) |>
  ## scale the variables
  mutate(y = ifelse(y > 4, 1, 0)) |>
  mutate(y_std = (y - mean(y)) / sd(y)) |>
  mutate(t = scales::rescale(t)) |>
  ## select and wrap up organization
  select(pid, t, y, y_std) |> arrange(pid, t)

# --- reference
d_slopes <- data |>
  nest(.by = "pid") |>
  mutate(m  = map(
    .x = data, 
    .f = ~ lm(y_std ~ t, data = .))) |>
  mutate(lm = map_dbl(
    .x = m, 
    .f = ~ round(coef(.)[[2]], 3)))  |>
  select(pid, lm)

# PART 2: ILLUSTRATION ------------------------------------------------------- #

set.seed(112358)
d_example <- buildDGP(
  n = data |> distinct(pid) |> nrow(),
  t = data |> distinct(t)   |> nrow(),
  ## change
  rate = 0.25, strength = 2,
  ## balance issues
  balance_dir = 0.5,
  balance_res = data |> pull(y) |> mean(),
  ## reliability
  reliable = 0.9,
  ## data exports
  export = FALSE, patterns = FALSE, slopes = TRUE)

# --- plot 1
example_p1 <- d_slopes |>
  ggplot(aes(x = lm, y = after_stat(ndensity))) +
  geom_histogram(
    bins = 15,
    color = "black",
    fill = my_oka[1],
    alpha = 0.8
  ) +
  theme_ipsum_rc(grid = "XxY", axis_title_size = 12) +
  labs(subtitle = "Observed Slopes", 
       x = "", 
       y = "Relative Frequency") +
  scale_x_continuous(
    limits = c(-3, 3),
    oob = scales::oob_keep,
    breaks = scales::pretty_breaks()
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks())

# --- plot 2
example_p2 <- d_example |>
  ggplot(aes(x = estimate, y = after_stat(ndensity))) +
  geom_histogram(
    bins = 15,
    color = "black",
    fill = my_oka[3],
    alpha = 0.8
  ) +
  theme_ipsum_rc(grid = "XxY", axis_title_size = 12) +
  labs(subtitle = "Simulated Slopes", 
       x = "Individual-Level Slope Coefficients", 
       y = "") +
  scale_x_continuous(
    limits = c(-3, 3),
    oob = scales::oob_keep,
    breaks = scales::pretty_breaks()
  ) +
  scale_y_continuous(breaks = scales::pretty_breaks())

# --- plot 3
example_p3 <-
  tibble(observed = d_slopes$lm, simulate = d_example$estimate) |> 
  ggplot() +
  stat_ecdf(aes(x = observed), geom = "step", col = my_oka[1]) +
  stat_ecdf(aes(x = simulate), geom = "step", col = my_oka[3]) +
  theme_ipsum_rc(grid = "XxY", axis_title_size = 12) +
  labs(subtitle = "Empirical CDF",
       x = "Individual-Level Slope Coefficients", 
       y = "ECDF") +
  scale_x_continuous(breaks = scales::pretty_breaks()) +
  scale_y_continuous(breaks = scales::pretty_breaks())

# --- draw trajectories
png(
  "figures/immigration01.png",
  w = 12,
  h =  7,
  units = "in",
  res = 500
)
(example_p1 / example_p2) | example_p3
dev.off()

# PART 3: FAKE-DATA SIMULATIONS ---------------------------------------------- #

# --- perform simulations
d <- gridsearch::gridSearch(
  ## datafile
  data = data,
  ## vars
  yname = 'y', tname = 't', pname = 'pid',
  ## patterns
  pattern = "slopes", steps = 1,
  ## reliability
  reliability = .90, ## assumption based on ICC value
  workers = 10, seed = 112358)

saveRDS(d, "./data/output/immigration.rds")

# --- plot the distributions
png(
  "figures/immigration02.png",
  w = 8,
  h = 6,
  units = "in",
  res = 500
)
gridPlot(
  data = d,
  cut1 = 0.01,
  cut2 = 0.05,
  cut3 = 0.10
)
dev.off()

# ---------------------------------------------------------------------------- #
