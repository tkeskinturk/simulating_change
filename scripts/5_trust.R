
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
  purrr,
  furrr,
  srvyr,
  hrbrthemes,
  gridsearch,
  tidyverse)
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
  haven::read_dta("./data/britishelection.dta") |>
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
    y = predicted,
    ymin = conf.low,
    ymax = conf.high
  )) +
  ## statistics
  geom_line(linetype = "dashed", linewidth = 0.2) +
  geom_pointrange(size = 0.5) +
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

# PART 2: FAKE-DATA SIMULATIONS ---------------------------------------------- #

# --- perform simulations
d <- gridsearch::gridSearch(
  ## datafile
  data = data,
  ## vars
  yname = 'y', tname = 't', pname = 'pid',
  ## patterns
  pattern = "contingency", steps = 30,
  ## reliability
  reliability = 0.88, ## assumption based on ICC value
  workers = 10, seed = 112358)

saveRDS(d, "./data/output/trust.rds")

# --- plot the distributions
png(
  "figures/trust02.png",
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
