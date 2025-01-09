
# - promises and pitfalls of panel change ------------------------------------ #
# - case 1 # grass ----------------------------------------------------------- #

### note: this script
###       (a) takes `grass` as an examplary GSS item,
###       (b) simulates potential DGPs approximating the observed parameters,
###       (c) shows the relative fit of the DGPs to observed data.

rm(list = ls()) # clean-up

# install.packages("devtools")
# devtools::install_github("tkeskinturk/gridsearch")

pacman::p_load(
  srvyr,
  gssr,
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

# PART 1: GSS DATA ----------------------------------------------------------- #

data("gss_panel06_long")
data("gss_panel08_long")
data("gss_panel10_long")

# --- get data
gss <-
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
  )

# --- clean-up
gss <- gss |>
  mutate(grass = as.character(grass)) |>
  mutate(grass = case_when(grass == "Legal" ~ 1,
                           grass == "Not Legal" ~ 0,
                           .default = NA_real_)) |>
  drop_na() |>
  mutate(n = n(), .by = pid) |> 
  filter(n == 3) |> 
  select(-n)

# --- draw trajectories
png(
  "figures/grass01.png",
  w = 7,
  h = 5,
  units = "in",
  res = 500
)
gss |>
  ## weight
  as_survey(weights = c(wtpannr123)) |>
  ## summarize
  group_by(year, wave, run) |>
  srvyr::summarize(p = survey_mean(grass)) |>
  ungroup() |>
  ## draw
  ggplot(aes(x = year, y = p, col = run)) +
  ## statistics
  geom_line(aes(group = run), linetype = "dashed", linewidth = 0.2) +
  geom_point(size = 2) +
  ## labels
  labs(x = "Year", y = "% Support", 
       col = "GSS Panel Cohort") +
  scale_y_continuous(
    limits = c(0.35, 0.61),
    labels = scales::label_percent()) +
  ## structure
  theme_ipsum_rc(grid = "XY", 
                 axis_title_size = 12) +
  theme(legend.position = "top") +
  ## colors
  scale_color_manual(values = c(my_oka[1:3]))
dev.off()

# PART 2: ABC SIMULATIONS ---------------------------------------------------- #

set.seed(112358) # for replicability

# --- grid search
d <-
  gridsearch::gridSearch(
    ## datafile
    data = gss,
    rel_min = .6,
    rel_max =  1,
    ## vars
    yname = 'grass',
    tname = 'wave',
    pname = 'pid',
    ## technicalities
    n_samples = 1e5
  )

saveRDS(d, "./data/output/grass.rds")

# PART 3: ABC DISTRIBUTIONS -------------------------------------------------- #

# --- plot the parameters
png(
  "figures/grass02.png",
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
  "figures/grass03.png",
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
