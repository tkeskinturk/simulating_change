
# - promises and pitfalls of panel change ------------------------------------ #
# - case study 2 # abortion -------------------------------------------------- #

### note: this script
###       (a) takes `abany` as an examplary GSS item,
###       (b) simulates potential DGPs approximating the observed parameters,
###       (c) shows the relative fit of the DGPs to observed data.

rm(list = ls()) # clean-up
pacman::p_load(
  gssr,
  ggh4x,
  purrr,
  furrr,
  hrbrthemes,
  srvyr,
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

# PART 1: GSS DATA ----------------------------------------------------------- #

data("gss_panel06_long")
data("gss_panel08_long")
data("gss_panel10_long")

# --- get data
gss <-
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

# --- draw trajectories
png(
  "figures/abany01.png",
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
  srvyr::summarize(p = survey_mean(abany)) |>
  ungroup() |>
  ## draw
  ggplot(aes(x = year, y = p, col = run)) +
  ## statistics
  geom_line(linetype = "dashed", linewidth = 0.2) +
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

# --- reference
gss_contingency <- gss |>
  select(pid, wave, abany) |>
  pivot_wider(
    names_from = "wave", values_from = "abany") |>
  janitor::clean_names() |>
  unite("patterns", starts_with("x"), sep = "") |>
  summarize(reference = n(), 
            .by = "patterns") |>
  arrange(patterns)

# PART 2: SIMS --------------------------------------------------------------- #

# --- build-up the parameter space
gss_sim <-
  expand_grid(
    p_n =
      gss |> distinct(pid) |> nrow(),
    p_t =
      3,
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
    p_res = gss |> pull(abany) |> mean(),
    p_rel = 0.86, # `abany` from Hout & Hastings
    sims = 1:1000
  )

# --- simulate panel datasets from the parameters
plan(multisession)
set.seed(112358)
gss_sim <- gss_sim |>
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
            export = FALSE, patterns = TRUE, slopes = FALSE
          )
        ),
        ## for reproducibility
        .options = furrr_options(seed = TRUE),
        .progress = TRUE
      )
  )
saveRDS(gss_sim, "./data/output/abany.rds")

# PART 3: MAPPING ------------------------------------------------------------ #

# --- perform the mapping
gss_sum <- gss_sim |>
  unnest(data) |>
  left_join(gss_contingency, by = "patterns") |>
  mutate(deviation = abs(sim_counts - reference)) |>
  summarize(
    deviation_sum = sum(deviation),
    .by = c("p_rate", "p_strength", "p_dir", "sims")
  ) |>
  mutate(deviation_sum = deviation_sum / gss |> distinct(pid) |> nrow()) |>
  summarize(error = mean(deviation_sum),
            .by = c("p_rate", "p_strength", "p_dir"))

# --- clean-up the dataframe
tag <- which(gss_sum$error == min(gss_sum$error))
gss_sum <- gss_sum |>
  ## cleanup
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
  mutate(
    group = case_when(
      error <= 0.05 & error <= quantile(error, 0.05) ~
        "Group 1",
      error <= 0.10 & error <= quantile(error, 0.10) ~
        "Group 2",
      error <= 0.20 & error <= quantile(error, 0.20) ~
        "Group 3",
      TRUE ~ "Group 4"
    )
  ) |>
  mutate(group = ifelse(row_number() == tag, "Group 5", group)) |>
  mutate(group = factor(
    group,
    levels = c("Group 1", "Group 2", "Group 3", "Group 4", "Group 5")
  ))

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
gss_plot <- gss_sum |> 
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
  "figures/abany02.png",
  w = 8,
  h = 6,
  units = "in",
  res = 500
)
gss_plot + ggh4x::facet_manual(~p_dir, design = design)
dev.off()

# ---------------------------------------------------------------------------- #
