
# - promises and pitfalls of panel change ------------------------------------ #
# - script 1 # parameter space for `gridSearch` ------------------------------ #

### note: this script
###       (a) generates simulated DGPs with varying parameters,
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
  broom.helpers)
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
    p_t = c(3, 5, 7),
    ## the rate of change in the population
    p_rate = c(.1, .25, .5),
    ## the direction of change (halfsies versus secular)
    p_dir = c(.5, 1),
    ## the fixed level of marginal distributions
    p_res = c(.5, .75),
    ## the strength of change in latent variable
    p_strength = c(.5, 1),
    ## the reliability score of the outcome measurement
    p_reliable = .9,
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

d <- d |>
  mutate(
    grid = map(.x = data,
               .f = ~gridSearch(data = .,
                                yname = "y_obs",
                                tname = "t",
                                pname = "pid",
                                pattern = "contingency",
                                steps = 1,
                                reliability = 0.9,
                                workers = 10,
                                seed = 11235)
    ))
saveRDS(d, file = "./data/simulations/simulations.rds")

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
  mutate(rank = future_pmap_dbl(
    .l = list(grid, p_rate, p_dir, p_strength),
    .f = ~ . |> 
      arrange(error) |>
      mutate(rank = row_number()) |> 
      filter(rate == ..2 & direction == ..3 & strength == ..4) |>
      pull(rank)
  )) |>
  ## clean-up
  mutate(
    time = factor(
      p_t,
      levels = c(3, 5, 7),
      labels = c("3 Waves", "5 Waves", "7 Waves")
    ),
    rate = factor(
      p_rate,
      levels = c(0.1, 0.25, 0.5),
      labels = c("10% Rate", "25% Rate", "50% Rate")
    ),
    direction = p_dir,
    marginals = factor(
      p_res,
      levels = c(0.5, 0.75),
      labels = c("50% Marginals", "75% Marginals")
    ),
    strength = factor(
      p_strength,
      levels = c(0.5, 1),
      labels = c("0.5 SD", "1 SD")
    )
  ) |>
  select(time, rate, strength, direction, marginals, sims, rank)

# PART 4: EXPLORE ------------------------------------------------------------ #

## pretty variables
variable_lookup <- tribble(
  ~ variable,
  ~ variable_nice,
  "time",
  "Time",
  "rate",
  "Rate of Change",
  "strength",
  "Strength",
  "direction",
  "Direction",
  "marginals",
  "Marginals"
) |>
  mutate(variable_nice = fct_inorder(variable_nice))

## calculate the coefficients
coefs <- lm(rank ~
              time +
              rate +
              direction +
              marginals +
              strength, data = d) |> 
  ## some andrew heiss magick!
  tidy_and_attach() |>
  tidy_add_reference_rows() |>
  tidy_add_estimate_to_reference_rows() |>
  filter(term != "(Intercept)") |>
  mutate(term_nice = str_remove(term, variable)) |>
  left_join(variable_lookup, by = join_by(variable)) |>
  mutate(across(c(term_nice, variable_nice), 
                ~ fct_inorder(.))) |>
  mutate(across(c("conf.low", "conf.high"),
                ~ ifelse(is.na(.) == TRUE, 0, .)))

## plot coefficients
png(
  "./figures/simulation01.png",
  w = 7,
  h = 7,
  units = "in",
  res = 500
)
ggplot(coefs, aes(x = estimate, 
                  y = fct_rev(term_nice))) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             col = "gray25") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high)) +
  labs(x = "Change in Rank", y = NULL) +
  facet_wrap(~ variable_nice, scales = "free_y", ncol = 1)
dev.off()

# ---------------------------------------------------------------------------- #
