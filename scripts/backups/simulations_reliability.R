
# - promises and pitfalls of panel change ------------------------------------ #
# - script 2 # reliability space for `gridSearch` ---------------------------- #

### note: this script
###       (a) generates simulated DGPs with varying reliability scores,
###       (b) applies the grid search algorithm to recover error distributions,
###       (c) evaluates whether the algorithm can recover "true" DGP.

rm(list = ls()) # clean-up

# install.packages("devtools")
# devtools::install_github("tkeskinturk/gridsearch")

pacman::p_load(
  purrr,
  furrr,
  hrbrthemes,
  gridsearch,
  ggforce,
  tidyverse,
  marginaleffects)
theme_set(theme_ipsum_rc())

# --- helpers
my_oka <-
  c("#eebd64",
    "#8b3b36", 
    "#a2b4b7")
`%nin%` = Negate(`%in%`)

# PART 1: DATA --------------------------------------------------------------- #

set.seed(11235)

# --- conditions
d <-
  expand_grid(
    ## 500 individuals
    p_n = 500,
    ## varying panel waves (restrict to contingency)
    p_t = 3,
    ## the rate of change in the population
    p_rate = c(.1, .25, .5),
    ## the direction of change (halfsies versus secular)
    p_dir = c(.5, 1),
    ## the fixed level of marginal distributions
    p_res = .5,
    ## the strength of change in latent variable
    p_strength = c(.5, 1),
    ## the reliability score of the outcome measurement
    p_reliable = c(.7, .9),
    ## the number of simulated datasets
    sims = 1:30)

# --- dataframes
d <- d |>
  dplyr::mutate(
    data =
      furrr::future_pmap(
        ## mapping list
        .l = list(p_n,
                  p_t,
                  p_rate,
                  p_dir,
                  p_res,
                  p_strength,
                  p_reliable),
        ## refer to list based on index
        .f = purrr::possibly(
          ~ buildDGP(
            # varying parameters
            n = ..1,
            t = ..2,
            rate = ..3,
            balance_dir = ..4,
            balance_res = ..5,
            strength = ..6,
            reliable = ..7,
            export = TRUE,
            patterns = FALSE, 
            slopes = FALSE
          )
        ),
        ## for reproducibility
        .options = furrr::furrr_options(seed = TRUE),
        .progress = TRUE
      )
  )

# PART 2: GRID SEARCH -------------------------------------------------------- #

# .70 reliability
d <- d |>
  mutate(
    grid_70 = map(
      .x = data,
      .f = ~ gridSearch(
        data = .,
        yname = "y_obs",
        tname = "t",
        pname = "pid",
        pattern = "contingency",
        steps = 1,
        reliability = 0.7,
        workers = 10,
        seed = 11235
      )
    ))

# .90 reliability
d <- d |>
  mutate(
    grid_90 = map(
      .x = data,
      .f = ~ gridSearch(
        data = .,
        yname = "y_obs",
        tname = "t",
        pname = "pid",
        pattern = "contingency",
        steps = 1,
        reliability = 0.9,
        workers = 10,
        seed = 11235
      )
    ))

saveRDS(d, file = "./data/simulations/reliability.rds")

# PART 3: ORGANIZATION ------------------------------------------------------- #

d <- d |>
  mutate(p_dir = factor(
    p_dir,
    levels = c(0, 0.25, 0.5, 0.75, 1),
    labels = c(
      "100% Down",
      "75% Down-25% Up",
      "50% Down-50% Up",
      "25% Down-75% Up",
      "100% Up"
    )
  )) |>
  ## DGP rank
  mutate(rank_70 = future_pmap_dbl(
    .l = list(grid_70, p_rate, p_dir, p_strength),
    .f = ~ . |>
      arrange(error) |>
      mutate(rank = row_number()) |>
      filter(rate == ..2 & direction == ..3 & strength == ..4) |>
      pull(rank)
  )) |>
  mutate(rank_90 = future_pmap_dbl(
    .l = list(grid_90, p_rate, p_dir, p_strength),
    .f = ~ . |>
      arrange(error) |>
      mutate(rank = row_number()) |>
      filter(rate == ..2 & direction == ..3 & strength == ..4) |>
      pull(rank)
  )) |>
  ## clean-up
  mutate(
    rate = factor(
      p_rate,
      levels = c(0.1, 0.25, 0.5),
      labels = c("10% Rate", "25% Rate", "50% Rate")
    ),
    direction = p_dir,
    strength = factor(
      p_strength,
      levels = c(0.5, 1),
      labels = c("0.5 SD", "1 SD")
    ),
    reliability = factor(
      p_reliable,
      levels = c(0.7, 0.9),
      labels = c("70%", "90%")
    )
  ) |>
  select(rate, strength, direction, reliability, sims, rank_70, rank_90) |>
  pivot_longer(
    cols = c(rank_70, rank_90),
    names_to = "grid",
    values_to = "rank"
  ) |>
  mutate(grid = factor(
    grid,
    levels = c("rank_70", "rank_90"),
    labels = c("70%", "90%")
  ))

# PART 4: EXPLORE ------------------------------------------------------------ #

## pretty variables
variable_lookup <- tribble(
  ~ variable,
  ~ variable_nice,
  "rate",
  "Rate of Change",
  "strength",
  "Strength",
  "direction",
  "Direction",
  "reliability_set",
  "Reliability"
) |>
  mutate(variable_nice = fct_inorder(variable_nice))

## cross-classify
d <- d |>
  mutate(reliability_set = paste(reliability, grid, sep = " - ")) |> 
  mutate(reliability_set = factor(
    reliability_set,
    levels = c("70% - 70%",
               "70% - 90%",
               "90% - 70%",
               "90% - 90%")
  ))

## calculate the coefficients
coefs <- lm(rank ~
              reliability_set +
              rate +
              direction +
              strength, data = d)

## plot coefficients
png(
  "./figures/simulation02.png",
  w = 8,
  h = 4,
  units = "in",
  res = 500
)
plot_predictions(
  model = coefs,
  condition = "reliability_set",
  vcov = sandwich::vcovCL(coefs, cluster = d$sims)
) +
  labs(x = "Condition", y = "Rank") +
  scale_y_continuous(limits = c(1200, 4800),
                     breaks = scales::pretty_breaks(n = 3))
dev.off()

# ---------------------------------------------------------------------------- #
