
# - promises and pitfalls of panel change ------------------------------------ #
# - case 3 # trust ----------------------------------------------------------- #

### note: this script
###       (a) takes `trust` as an examplary BES item,
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

# PART 1: BES DATA ----------------------------------------------------------- #

# --- get data
data <-
  haven::read_dta("./data/source_files/bes.dta") |>
  select(
    pid = id,
    trust_1 = genTrustW17,
    trust_2 = genTrustW20,
    trust_3 = genTrustW21,
    trust_4 = genTrustW23,
    trust_5 = genTrustW27
  ) |>
  haven::zap_labels() |>
  mutate_at(
    .vars = vars(matches("trust")),
    .funs = ~ case_when(. == 1 ~ 1, . == 2 ~ 0, TRUE ~ NA_real_)
  ) |>
  drop_na() |>
  pivot_longer(cols = -pid,
               names_to = "var",
               values_to = "y") |>
  separate(var, into = c("var", "t")) |>
  mutate(t = as.integer(t)) |> select(-var) |> 
  mutate(
    year =
      case_when(
        t == 1 ~ 2019,
        t == 2 ~ 2020,
        t == 3 ~ 2021,
        t == 4 ~ 2022,
        t == 5 ~ 2024
      )
  )

# --- draw trajectories
png(
  "figures/trust01.png",
  w = 7,
  h = 5,
  units = "in",
  res = 500
)
glm(y ~ factor(year), data = data, family = "binomial") |>
  ggeffects::ggemmeans(terms = "year") |>
  as_tibble() |>
  ## draw
  ggplot(aes(
    x = x,
    y = predicted
  )) +
  ## statistics
  geom_line(linetype = "dashed", linewidth = 0.2) +
  geom_point(size = 2) +
  ## labels
  labs(x = "Time", y = "% Trust") +
  scale_y_continuous(
    limits = c(0.30, 0.60),
    labels = scales::label_percent()) +
  ## structure
  theme_ipsum_rc(grid = "XY", 
                 axis_title_size = 12) +
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

saveRDS(d, "./data/output/trust.rds")

# PART 3: ABC DISTRIBUTIONS -------------------------------------------------- #

# --- plot the parameters
png(
  "figures/trust02.png",
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
  "figures/trust03.png",
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
