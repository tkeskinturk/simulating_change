
# - promises and pitfalls of panel change ------------------------------------ #
# - case 4 # immigration ----------------------------------------------------- #

### note: this script
###       (a) takes `imm` as an examplary panel data item,
###       (b) simulates potential DGPs approximating the observed parameters,
###       (c) shows the relative fit of the DGPs to observed data.

rm(list = ls()) # clean-up

# install.packages("devtools")
# devtools::install_github("tkeskinturk/gridsearch")

pacman::p_load(
  srvyr,
  gridsearch,
  tidyverse,
  hrbrthemes,
  GGally)
theme_set(theme_ipsum_rc())

# --- helpers
my_oka <-
  c("#eebd64",
    "#8b3b36", 
    "#a2b4b7")
`%nin%` = Negate(`%in%`)

# PART 1: 26-WAVE DATA ------------------------------------------------------- #

# --- get data
data <- 
  read_csv("./data/source_files/yllanon.csv") |> 
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
  ## select and wrap up organization
  select(pid, t, y) |> arrange(pid, t)

# --- draw trajectories
png(
  "figures/immigration01.png",
  w = 7,
  h = 5,
  units = "in",
  res = 500
)
data |>
  group_by(t) |>
  summarize(y = mean(y)) |>
  ## draw
  ggplot(aes(x = t, y = y)) +
  ## statistics
  geom_line(linetype = "dashed", linewidth = 0.2) +
  geom_point(size = 2) +
  ## labels
  labs(x = "Time", y = "Immigration Spending Preference") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 17)) +
  scale_y_continuous(
    limits = c(0.35, 0.55),
    labels = scales::label_percent()) +
  ## structure
  theme_ipsum_rc(grid = "XY", axis_title_size = 12) +
  theme(legend.position = "top")
dev.off()

# PART 2: ABC SIMULATIONS ---------------------------------------------------- #

set.seed(112358) # for replicability

# --- grid search
d <-
  gridsearch::gridSearch(
    ## datafile
    data = data,
    rel_min = .6,
    rel_max =  1,
    ## vars
    yname = 'y',
    tname = 't',
    pname = 'pid',
    ## technicalities
    n_samples = 1e5
  )

saveRDS(d, "./data/output/immigration.rds")

# PART 3: ABC DISTRIBUTIONS -------------------------------------------------- #

# --- plot the parameters
png(
  "figures/immigration02.png",
  w = 10,
  h =  8,
  units = "in",
  res = 500
)
d |>
  filter(error < quantile(d$error, .02)) |>
  ggpairs(
    columns = 1:4,
    columnLabels = c("Strength", "Rate", "Balance", "Reliability"),
    lower = list(continuous = wrap(
      "points", alpha = 0.5, size = .75
    ))
  )
dev.off()

# --- plot the distributions
png(
  "figures/immigration03.png",
  w = 8,
  h = 6,
  units = "in",
  res = 500
)
d |>
  mutate(balance = paste("Balance", 
                         "\U2248", 
                         round(bal_sample, 1))) |>
  mutate(balance = as.factor(balance)) |>
  filter(error < quantile(d$error, .02)) |> # keep the best 2 percent
  ggplot(aes(x = pc_sample, y = ic_sample)) +
  geom_point(alpha = .5) +
  facet_wrap(~ balance, drop = FALSE) +
  labs(x = "Rate of Change", y = "Strength of Change") +
  scale_x_continuous(limits = c(0, 1), labels = scales::label_percent()) +
  scale_y_continuous(limits = c(0, 2), labels = scales::pretty_breaks())
dev.off()

# ---------------------------------------------------------------------------- #
