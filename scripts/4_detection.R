
# - promises and pitfalls of panel change ------------------------------------ #
# - change detection --------------------------------------------------------- #

### note: this script
###       (a) takes `grass`, `abany`, and `imm` as variable examples,
###       (b) simulates potential DGPs approximating the observed parameters,
###       (c) assesses whether we have leverage to detect individuals.

rm(list = ls()) # clean-up
pacman::p_load(
  gssr,
  ggh4x,
  purrr,
  furrr,
  hrbrthemes,
  srvyr,
  patchwork,
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

# --- prep
plan(multisession)
set.seed(112358)

# PART 1: DATA --------------------------------------------------------------- #

data("gss_panel06_long")
data("gss_panel08_long")
data("gss_panel10_long")

# --- grass
grass <-
  ## gather
  bind_rows(
    gss_panel06_long |>
      select(pid = firstid, year, wave, wtpannr123, grass) |>
      haven::zap_labels() |>
      mutate(pid = as.integer(paste0(
        2006, as.integer(pid)
      ))) |> 
      mutate(run = "2006-2010"),
    gss_panel08_long |>
      select(pid = firstid, year, wave, wtpannr123, grass) |>
      haven::zap_labels() |>
      mutate(pid = as.integer(paste0(
        2008, as.integer(pid)
      ))) |> 
      mutate(run = "2008-2012"),
    gss_panel10_long |>
      select(pid = firstid, year, wave, wtpannr123, grass) |>
      haven::zap_labels() |>
      mutate(pid = as.integer(paste0(
        2010, as.integer(pid)
      ))) |> 
      mutate(run = "2010-2014")
  ) |>
  ## organize
  drop_na() |>
  mutate(n = n(), .by = pid) |> filter(n == 3) |> select(-n) |>
  mutate(grass = ifelse(grass == 2, 0, 1))

# --- abany
abany <-
  ## gather
  bind_rows(
    gss_panel06_long |>
      select(pid = firstid, year, wave, wtpannr123, abany) |>
      haven::zap_labels() |>
      mutate(pid = as.integer(paste0(
        2006, as.integer(pid)
      ))) |> 
      mutate(run = "2006-2010"),
    gss_panel08_long |>
      select(pid = firstid, year, wave, wtpannr123, abany) |>
      haven::zap_labels() |>
      mutate(pid = as.integer(paste0(
        2008, as.integer(pid)
      ))) |> 
      mutate(run = "2008-2012"),
    gss_panel10_long |>
      select(pid = firstid, year, wave, wtpannr123, abany) |>
      haven::zap_labels() |>
      mutate(pid = as.integer(paste0(
        2010, as.integer(pid)
      ))) |> 
      mutate(run = "2010-2014")
  ) |>
  ## organize
  drop_na() |>
  mutate(n = n(), .by = pid) |> filter(n == 3) |> select(-n) |>
  mutate(abany = ifelse(abany == 2, 0, 1))

# --- immigration
immig <- 
  read_csv("./data/yllanon.csv") |> 
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

# PART 2: GRASS -------------------------------------------------------------- #

# --- parameter space
grass_sims_best <-
  expand_grid(
    p_n =
      grass |> distinct(pid) |> nrow(),
    p_t =
      3,
    p_rate = 0.25,
    p_strength = 0.9,
    p_dir = 1,
    p_res = grass |> pull(grass) |> mean(),
    p_rel = 0.91, # `grass` from Hout & Hastings
    sims = 1:1000
  )

# --- simulate panel datasets from the parameters
grass_sims_best <- grass_sims_best |>
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
            export = TRUE, patterns = FALSE, slopes = TRUE
          )
        ),
        ## for reproducibility
        .options = furrr_options(seed = TRUE),
        .progress = TRUE
      )
  )

# --- confusion matrices
grass_sims_best <- grass_sims_best |>
  unnest_wider(col = data) |>
  # step 1: classify
  mutate(correspondence = map2(
    .x = data,
    .y = slopes,
    .f =
      ~ left_join(.x |>
                    distinct(pid, ever_change), .y, by = "pid") |>
      mutate(predicted =
               ifelse(
                 abs(estimate - 0) < abs(estimate - 0.9), 0, 1
               )) |>
      select(pid, ever_change, predicted)
  )) |>
  # step 2: confuse
  mutate(confusion = map(
    .x = correspondence,
    .f = ~  tibble(
      ever_change =
        c(0, 0, 1, 1),
      predicted =
        c(0, 1, 0, 1)
    ) |>
      left_join(
        . |>
          count(ever_change, predicted),
        by = c("ever_change", "predicted")
      ) |>
      mutate(n = ifelse(is.na(n) == TRUE, 0, n))
  )) |>
  # step 3: organize
  mutate(confused =
           map(.x = confusion, .f = ~ confused_outcomes(.))) |>
  select(sims, confused) |> unnest_wider(confused)

# PART 3: ABANY -------------------------------------------------------------- #

# --- parameter space
abany_sims_best <-
  expand_grid(
    p_n =
      abany |> distinct(pid) |> nrow(),
    p_t =
      3,
    p_rate = 0.8,
    p_strength = 0.4,
    p_dir = 0.5,
    p_res = abany |> pull(abany) |> mean(),
    p_rel = 0.86, # `abany` from Hout & Hastings
    sims = 1:1000
  )

# --- simulate panel datasets from the parameters
abany_sims_best <- abany_sims_best |>
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
            export = TRUE, patterns = FALSE, slopes = TRUE
          )
        ),
        ## for reproducibility
        .options = furrr_options(seed = TRUE),
        .progress = TRUE
      )
  )

# --- confusion matrices
abany_sims_best <- abany_sims_best |>
  unnest_wider(col = data) |>
  # step 1: classify
  mutate(correspondence = map2(
    .x = data,
    .y = slopes,
    .f =
      ~ left_join(.x |>
                    distinct(pid, ever_change), .y, by = "pid") |>
      mutate(predicted =
               ifelse(
                 abs(estimate) < 0.4, 0, 1
               )) |>
      select(pid, ever_change, predicted)
  )) |>
  # step 2: confuse
  mutate(confusion = map(
    .x = correspondence,
    .f = ~  tibble(
      ever_change =
        c(0, 0, 1, 1),
      predicted =
        c(0, 1, 0, 1)
    ) |>
      left_join(
        . |>
          count(ever_change, predicted),
        by = c("ever_change", "predicted")
      ) |>
      mutate(n = ifelse(is.na(n) == TRUE, 0, n))
  )) |>
  # step 3: organize
  mutate(confused =
           map(.x = confusion, .f = ~ confused_outcomes(.))) |>
  select(sims, confused) |> unnest_wider(confused)

# PART 4: IMMIGRATION -------------------------------------------------------- #

# --- parameter space
immig_sims_best <-
  expand_grid(
    p_n =
      immig |> distinct(pid) |> nrow(),
    p_t =
      17,
    p_rate = 0.6,
    p_strength = 0.5,
    p_dir = 0,
    p_res = immig |> pull(y) |> mean(),
    p_rel = 0.9,
    sims = 1:1000
  )

# --- simulate panel datasets from the parameters
immig_sims_best <- immig_sims_best |>
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
            export = TRUE, patterns = FALSE, slopes = TRUE
          )
        ),
        ## for reproducibility
        .options = furrr_options(seed = TRUE),
        .progress = TRUE
      )
  )

# --- confusion matrices
immig_sims_best <- immig_sims_best |>
  unnest_wider(col = data) |>
  # step 1: classify
  mutate(correspondence = map2(
    .x = data,
    .y = slopes,
    .f =
      ~ left_join(.x |>
                    distinct(pid, ever_change), .y, by = "pid") |>
      mutate(predicted =
               ifelse(
                 abs(estimate - 0) < abs(estimate - (-0.5)), 0, 1
               )) |>
      select(pid, ever_change, predicted)
  )) |>
  # step 2: confuse
  mutate(confusion = map(
    .x = correspondence,
    .f = ~  tibble(
      ever_change =
        c(0, 0, 1, 1),
      predicted =
        c(0, 1, 0, 1)
    ) |>
      left_join(
        . |>
          count(ever_change, predicted),
        by = c("ever_change", "predicted")
      ) |>
      mutate(n = ifelse(is.na(n) == TRUE, 0, n))
  )) |>
  # step 3: organize
  mutate(confused =
           map(.x = confusion, .f = ~ confused_outcomes(.))) |>
  select(sims, confused) |> unnest_wider(confused)

# PART 5: DRAW --------------------------------------------------------------- #

d <- bind_rows(
  abany_sims_best |> mutate(item = "Abortion"),
  grass_sims_best |> mutate(item = "Grass"),
  immig_sims_best |> mutate(item = "Immigration")) |> 
  mutate(item = factor(item,
                       levels = c("Grass",
                                  "Abortion",
                                  "Immigration"))) |> 
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
  scale_y_continuous(limits = c(0.1, 0.6),
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
  scale_y_continuous(limits = c(0.5, 1),
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
