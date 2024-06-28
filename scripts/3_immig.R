
# - promises and pitfalls of panel change ------------------------------------ #
# - case study 3 # immigration ----------------------------------------------- #

### note: this script
###       (a) takes `imm` as an examplary panel data item,
###       (b) simulates potential DGP distributions using fixed-effects models,
###       (c) shows the relative fit of the DGPs to observed data.

rm(list = ls()) # clean-up
pacman::p_load(
  ggh4x,
  purrr,
  furrr,
  patchwork,
  hrbrthemes,
  tidyverse)
theme_set(theme_ipsum_rc())

# --- helpers
my_oka <-
  c("#eebd64",
    "#8b3b36", 
    "#a2b4b7")
`%nin%` = Negate(`%in%`)

# --- fire power
source("./scripts/helper_functions/build_DGP.R")
source("./scripts/helper_functions/classification_metrics.R")

# PART 1: 26-WAVE DATA ------------------------------------------------------- #

# --- read in
d <- read_csv("./data/yllanon.csv")

# --- cleaning up
d <- d |>
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
d_slopes <- d |>
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
d_example <- generate_data(
  n = d |> distinct(pid) |> nrow(),
  t = d |> distinct(t)   |> nrow(),
  ## change
  rate = 0.25, strength = 2,
  ## balance issues
  balance_dir = 0.5,
  balance_res = d |> pull(y) |> mean(),
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

# PART 3: SIMS --------------------------------------------------------------- #

# --- build-up the parameter space
d_sim <-
  expand_grid(
    p_n =
      d |> distinct(pid) |> nrow(),
    p_t =
      d |> distinct(t)   |> nrow(),
    p_rate =     seq(from = 0, 
                     to   = 1, 
                     length.out = 21),
    p_strength = seq(from = 0.1, 
                     to   = 2, 
                     length.out = 20),
    p_dir = c(0, 
              0.25, 
              0.50, 
              0.75, 
              1),
    p_res = d |> pull(y) |> mean(),
    p_rel = 0.9,
    sims = 1:1000
  )

# --- simulate panel datasets from the parameters
plan(multisession)
set.seed(112358)
d_sim <- d_sim |>
  mutate(
    data =
      furrr::future_pmap(
        ## mapping list
        .l = list(
          p_n,
          p_t,
          p_strength,
          p_rate,
          p_dir,
          p_res,
          p_rel),
        ## refer to list based on index
        .f = purrr::possibly(
          ~ generate_data(
            # varying parameters
            n = ..1,
            t = ..2,
            strength = ..3,
            rate = ..4,
            balance_dir = ..5,
            balance_res = ..6,
            reliable = ..7,
            export = FALSE, patterns = FALSE, slopes = TRUE
          )
        ),
        ## for reproducibility
        .options = furrr_options(seed = TRUE),
        .progress = TRUE
      )
  )
saveRDS(d_sim, file = "./data/output/immigrant.rds")

# PART 4: MAPPING ------------------------------------------------------------ #

# --- perform the mapping
d_sum <- d_sim |>
  mutate(ks = future_map(
    .x = data,
    .f = ~ ks.test(
      x = .$estimate,
      y = d_slopes$lm,
      exact = TRUE
    )$statistic,
    .progress = TRUE
  )) |>
  unnest(cols = ks) |>
  summarize(ks_stat = mean(ks),
            .by = c("p_rate", "p_strength", "p_dir"))

# --- clean-up the dataframe
tag <- which(d_sum$ks_stat == min(d_sum$ks_stat))
d_sum <- d_sum |>
  ## cleanup
  mutate(p_dir = factor(p_dir,
                        levels = c(0, 0.25, 0.5, 0.75, 1),
                        labels = c("100% Down",
                                   "75% Down-25% Up",
                                   "50% Down-50% Up",
                                   "25% Down-75% Up",
                                   "100% Up"))) |>
  mutate(
    group = case_when(
      ks_stat <= 0.05 & ks_stat <= quantile(ks_stat, 0.05) ~
        "Group 1",
      ks_stat <= 0.10 & ks_stat <= quantile(ks_stat, 0.10) ~
        "Group 2",
      ks_stat <= 0.20 & ks_stat <= quantile(ks_stat, 0.20) ~
        "Group 3",
      TRUE ~ "Group 4")
    ) |>
  mutate(
    group = ifelse(row_number() == tag, "Group 5", group)
    ) |>
  mutate(
    group = factor(group,
                   levels = c("Group 1", 
                              "Group 2", 
                              "Group 3", 
                              "Group 4", 
                              "Group 5")))

# --- draw

## facet design
design <- 
  c(
  "
  AABBCC
  #DDEE#
  "
  )

## plot
d_plot <- d_sum |> 
  ggplot(aes(x = p_rate, y = p_strength)) +
  scale_x_continuous(limits = c(0, 1),
                     labels = scales::label_percent()) +
  scale_y_continuous(limits = c(0, 2),
                     labels = scales::pretty_breaks()) +
  labs(x = "Rate of Change", y = "Strength of Change") +
  theme_ipsum_rc(grid = "XY") +
  geom_point(aes(color = group), size = 2) +
  scale_color_manual(values = c("#000000",
                                "#bbbbbb",
                                "#eeeeee",
                                "transparent",
                                my_oka[1]),
                     drop = FALSE) +
  theme(legend.position = "none")

## design
png(
  "figures/immigration02.png",
  w = 8,
  h = 6,
  units = "in",
  res = 500
)
d_plot + ggh4x::facet_manual(~p_dir, design = design)
dev.off()

# ---------------------------------------------------------------------------- #
