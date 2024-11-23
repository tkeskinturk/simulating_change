
# - promises and pitfalls of panel change ------------------------------------ #
# - change detection --------------------------------------------------------- #

### note: this script
###       (a) takes `grass`, `abany`, and `imm` as variable examples,
###       (b) simulates potential DGPs approximating the observed parameters,
###       (c) assesses whether we have leverage to detect individuals.

rm(list = ls()) # clean-up

# devtools (install if necessary!)
# devtools::install_github("tkeskinturk/gridsearch")

pacman::p_load(
  purrr,
  furrr,
  hrbrthemes,
  gridsearch,
  tidyverse,
  patchwork)
theme_set(theme_ipsum_rc())

# --- helpers
my_oka <-
  c("#eebd64",
    "#8b3b36", 
    "#a2b4b7",
    "#4c6337")
`%nin%` = Negate(`%in%`)

# PART 1: DATA --------------------------------------------------------------- #

grass <- readRDS("./data/output/grass.rds")
abany <- readRDS("./data/output/abany.rds")
trust <- readRDS("./data/output/trust.rds")
immig <- readRDS("./data/output/immigration.rds")

grass <- grass |> arrange(error) |> slice(1)
abany <- abany |> arrange(error) |> slice(1)
trust <- trust |> arrange(error) |> slice(1)
immig <- immig |> arrange(error) |> slice(1)

# --- grass
grass_sim <-
  simulateChangers(
    n = grass |> pull(n),
    t = 3,
    rate = grass |> pull(rate),
    balance_dir = 1,
    balance_res = 0.47,
    strength = grass |> pull(strength),
    reliable = 0.91
  )

# --- abany
abany_sim <-
  simulateChangers(
    n = abany |> pull(n),
    t = 3,
    rate = abany |> pull(rate), ## for help
    balance_dir = 1,
    balance_res = 0.45,
    strength = abany |> pull(strength),
    reliable = 0.86
  )

# --- trust
trust_sim <-
  simulateChangers(
    n = trust |> pull(n),
    t = 5,
    rate = trust |> pull(rate),
    balance_dir = 0.5,
    balance_res = 0.43,
    strength = trust |> pull(strength),
    reliable = 0.88
  )

# --- imm
immig_sim <-
  simulateChangers(
    n = immig |> pull(n),
    t = 17,
    rate = immig |> pull(rate),
    balance_dir = 0,
    balance_res = 0.44,
    strength = immig |> pull(strength),
    reliable = 0.90
  )

# PART 2: DRAW --------------------------------------------------------------- #

d <- bind_rows(
  grass_sim |> mutate(item = "Grass"),
  abany_sim |> mutate(item = "Abortion*"),
  trust_sim |> mutate(item = "Trust"),
  immig_sim |> mutate(item = "Immigration")
) |>
  mutate(item = factor(item, levels = c(
    "Grass", "Abortion*", "Trust", "Immigration"
  ))) |> 
  pivot_longer(cols = -c(sims, item)) |> 
  filter(name != "kappa") |>
  mutate(name = factor(name,
                       levels = c("sensitivity",
                                  "specificity"),
                       labels = c("Sensitivity",
                                  "Specificity")))

## plot 1
p1 <- d |>
  filter(name == "Sensitivity") |>
  ggplot(aes(x = item, y = value, col = item)) +
  geom_point(shape = 95,
             size = 20,
             alpha = .1) +
  scale_color_manual(values = my_oka) +
  theme_ipsum_rc(grid = "XY") + theme(legend.position = "none") +
  labs(x = "", y = "Sensitivity", title = "Sensitivity") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = scales::pretty_breaks(n = 6))

## plot 2
p2 <- d |>
  filter(name == "Specificity") |>
  ggplot(aes(x = item, y = value, col = item)) +
  geom_point(shape = 95,
             size = 20,
             alpha = .1) +
  scale_color_manual(values = my_oka) +
  theme_ipsum_rc(grid = "XY") + theme(legend.position = "none") +
  labs(x = "Case Studies", y = "Specificity", title = "Specificity") +
  scale_y_continuous(limits = c(0, 1),
                     breaks = scales::pretty_breaks(n = 6))

png(
  "./figures/detection.png",
  w = 10,
  h =  8,
  units = "in",
  res = 500
)
p1 / p2
dev.off()

# ---------------------------------------------------------------------------- #
