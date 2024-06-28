
# - promises and pitfalls of panel change ------------------------------------ #
# - helper: build DGP -------------------------------------------------------- #

### note: this script 
###       (a) generates panel data from specified parameters values,
###       (b) estimates user-specified models of change,
###       (c) stores (1) patterns, (2) model results OR (3) data, if specified.

require(tidyverse)

# PART 1: FUNCTION ----------------------------------------------------------- #

generate_data <-
  
  function(
    
    # ----------------------------------------------------- #
    # FUNCTION ARGUMENTS AND WARNINGS                       #
    # ----------------------------------------------------- #
    
    ## the number of Ns
    n = 1000,
    ## the number of Ts
    t = 3,
    ## the rate of change in the population
    rate = 0.25,
    ## % balance in changes
    balance_dir = 0.5,
    ## % balance in `y` = 1
    balance_res = 0.5,
    ## strength of change in latent scores
    strength = 1,
    ## reliability for latent scores
    reliable = 0.8,
    ## return options for components
    export = TRUE, patterns = TRUE, slopes = TRUE) {
    
    if(n <= 0)
      return(stop(
        "Error: We need at least *some* people. Check your `n` call."))
    if(t <= 1)
      return(stop(
        "Error: We need at least 2 time periods to generate panel data. Try again."))
    if(rate < 0)
      return(stop(
        "Error: Rate of change must be 0 or a positive number."))
    if(balance_dir <  0 | balance_dir >  1)
      return(stop(
        "Error: Directionality balance must be a number between 0 and 1."))
    if(balance_res <= 0 | balance_res >= 1)
      return(stop(
        "Error: Outcome balance must be between 0 and 1 (not inclusive)."))
    if(strength < 0)
      return(stop(
        "Error: Change strength can't be negative (remember: it's 'strength')"))
    if(reliable < 0)
      return(stop(
        "Error: Reliability can't be negative. Try again."))
    if(export == FALSE & patterns == FALSE & slopes == FALSE)
      return(stop(
        "You said no data, no patterns, or no model information. What do you want?"))
    
    # ----------------------------------------------------- #
    # BUILDING BLOCKS                                       #
    # ----------------------------------------------------- #
    
    # --- part 1: generate actors
    if (rate == 0) {
      u <- 
        tibble(
          ## actors
          pid = c(1:n),
          ## central tendency
          u = rnorm(n = n, mean = 0, sd = 1)
        )
    } else {
      u <-
        tibble(
          ## actors
          pid = c(1:n),
          ## time at which they go through change
          change_time = rexp(n, rate = rate),
          ## central tendency
          u = rnorm(n = n, mean = 0, sd = 1)
        )
    }
    
    # --- part 2: generate the longitudinal data setup

    ## 2.1: build the window
    span_window <- seq(
      ## first t
      from = 0,
      ## final t                 
      to =   1,
      ## the number of ts
      length.out = t)
    
    ## 2.2: generate pid x time grid
    data <- expand_grid(u, t = span_window)
    
    ## 2.3: assign "change" status
    if (rate == 0) {
      data <- data |>
        mutate(change = 0) # placeholder for few clutter
    } else {
      data <- data |>
        mutate(change = ifelse(change_time <= span_window, 1, 0)) |>
        select(-change_time)
    }
    
    # ----------------------------------------------------- #
    # RESPONSE CONSTRUCTION                                 #
    # ----------------------------------------------------- #
    
    # --- part 3: generate the observed scores
    
    ## 3.1: add change scores
    data <- data |> 
      left_join(u |> select(pid) |>
                  mutate(upper = rbinom(
                    n = n, size = 1, prob = balance_dir
                  )) |>
                  mutate(upper = ifelse(upper == 1, 1, -1)),
                by = "pid"
      ) |> 
      mutate(u = u + strength * change * upper) |>
      select(-upper)
    
    ## 3.2: extract the new error variance
    new_error <- sd(data$u)
    
    ## 3.3: generate the realized y scores
    data <- data |>
      mutate(y =
               ## true scores
               (u * sqrt(reliable))
             +
               ## error
               rnorm(n = n, mean = 0, sd = new_error) 
             * 
               sqrt((1 - reliable)))
    
    # --- part 4: binarize response
    data <- data |>
      mutate(y_obs =
               ifelse(y <= quantile(y, prob = 1 - balance_res), 
                      0, 
                      1))
    
    # --- part 5: organize and spit out if necessary
    
    ## 5.1: organize the data
    data <- data |>
      mutate(ever_change = 
               ifelse(mean(change) > 0, 1, 0), 
             .by = "pid") |>
      dplyr::select(pid,
                    ever_change,
                    t, 
                    y_true = u, 
                    y_obs)
    
    ## 5.2: spit out the data
    if (export   == TRUE   & 
        patterns == FALSE  & 
        slopes   == FALSE) {
      return(data)
    }
    
    # ----------------------------------------------------- #
    # PATTERNS                                              #
    # ----------------------------------------------------- #
    
    data_patterns <- data |>
      select(pid, t, y_obs) |>
      pivot_wider(names_from = "t", values_from = "y_obs") |>
      janitor::clean_names() |>
      unite("patterns", starts_with("x"), sep = "") |>
      summarize(sim_counts = n(), .by = "patterns") |>
      arrange(patterns)
    
    if (export   == TRUE   & 
        patterns == TRUE   & 
        slopes   == FALSE) {
      output <- list(data = data, 
                     patterns = data_patterns)
      return(output)
    }
    
    if (export   == FALSE  & 
        patterns == TRUE   & 
        slopes   == FALSE) {
      return(data_patterns)
    }

    # ----------------------------------------------------- #
    # SLOPES                                                #
    # ----------------------------------------------------- #
    
    # standardize for `d` comparison
    data <- data |> 
      mutate(y_std = (y_obs - mean(y_obs)) / sd(y_obs))
    
    # linear regression
    ## model matrix
    model <- data |> nest(.by = "pid") |>
      mutate(m = map(.x = data, 
                     .f = ~ lm.fit(x = model.matrix(~ t, data = .), y = .$y_std)))
    ## extract estimates
    estimates_coef <-
      model |> mutate(estimate =
                        map_dbl(.x = m, 
                                .f = ~ round(coef(.)[[2]], 3))) |>
      select(pid, estimate)
    
    # --- part M3: export
    if (export   == FALSE  & 
        patterns == FALSE  & 
        slopes   == TRUE)  {
      return(estimates_coef)
    }
    
    if (export   == FALSE  & 
        patterns == TRUE   & 
        slopes   == TRUE)  {
      output <- list(patterns = data_patterns, 
                     slopes = estimates_coef)
      return(output)
    }
    
    if (export   == TRUE   & 
        patterns == TRUE   & 
        slopes   == TRUE)  {
      output <- list(data = data, 
                     patterns = data_patterns, 
                     slopes = estimates_coef)
      return(output)
    }
    
    if (export   == TRUE   & 
        patterns == FALSE  & 
        slopes   == TRUE)  {
      output <- list(data = data, 
                     slopes = estimates_coef)
      return(output)
    }

  }

# ---------------------------------------------------------------------------- #
