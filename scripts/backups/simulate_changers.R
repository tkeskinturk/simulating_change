
# - `gridsearch` ------------------------------------------------------------- #
# - function: simulate changers ---------------------------------------------- #

### note: this function
###       (a) takes a set of parameters from a specific DGP and simulates data,
###       (b) calculates accuracy scores for detecting changers.

# FUNCTION CALL -------------------------------------------------------------- #

#' Simulate Changers and Calculate Accuracy
#'
#' This function uses parameter values from a given DGP to calculate kappa, sensitivity, and specificity scores in simulated data frames.
#' @param n The number of individuals in the DGP.
#' @param t The number of time periods.
#' @param rate The percent of individuals changing across the panel.
#' @param balance_dir The direction of change, where 0 codes negative change, 1 codes positive change, and all values in-between codes the percentage of changers changing in the positive direction.
#' @param balance_res The marginal distribution of the outcome (effectively regulating the percent distribution of 0s and 1s).
#' @param strength The strength of change in the latent variable.
#' @param reliable The reliability score of the outcome measurement.
#' @param nrep The number of simulation runs.
#' @param seed Seed for reproducibility.
#' @param workers The number of workers for parallelization (note that parallelization is highly recommended for reasonable duration).
#'
#' @return A data frame.
#'
#' @details # Examples
#'
#' **Basic [simulateChangers()] call:**
#' ```{r, comment = "#>", collapse = TRUE}
#' set.seed(11235)
#' data <- simulateChangers(n = 100,
#'                          t = 3,
#'                          rate = 0.5,
#'                          balance_dir = 1,
#'                          balance_res = 0.5,
#'                          strength = 1,
#'                          reliable = 0.9,
#'                          nrep = 100,
#'                          seed = 11234,
#'                          workers = 10)
#' print(data)
#' ```
#'
#' @export
#'

simulateChangers <-

  function(

    # ----------------------------------------------------- #
    # FUNCTION ARGUMENTS AND WARNINGS                       #
    # ----------------------------------------------------- #

    n = 1000,
    t = 3,
    rate = 0.25,
    balance_dir = .5,
    balance_res = .5,
    strength = .5,
    reliable = 1,
    nrep = 1000, seed = 11235, workers = 8

  ) {

    if (n <= 0) stop("Error: We need at least some people. Check your `n` call.")
    if (t <= 1) stop("Error: We need at least 2 time periods to generate panel data. Try again.")
    if (rate < 0 | rate > 1) stop("Error: Rate of change must range between 0 and 1.")
    if (balance_dir <  0 | balance_dir >  1) stop("Error: Directionality balance must be a number between 0 and 1.")
    if (balance_res <= 0 | balance_res >= 1) stop("Error: Outcome balance must be between 0 and 1 (not inclusive).")
    if (strength < 0) stop("Error: Change strength can't be negative (remember: it's 'strength')")
    if (reliable < 0) stop("Error: Reliability can't be negative. Try again.")

    future::plan(future::multisession, workers = workers)
    set.seed(seed)

    # ----------------------------------------------------- #
    # BUILDING BLOCKS                                       #
    # ----------------------------------------------------- #

    # --- parameter space

    data <-
      tidyr::expand_grid(
        p_n = n,
        p_t = t,
        p_rate = rate,
        p_dir = balance_dir,
        p_res = balance_res,
        p_strength = strength,
        p_rel = reliable,
        sims = 1:nrep
      )

    # --- simulate panel datasets from the parameters
    data <- data |>
      dplyr::mutate(
        data =
          furrr::future_pmap(
            ## mapping list
            .l = list(.data$p_n,
                      .data$p_t,
                      .data$p_strength,
                      .data$p_rate,
                      .data$p_dir,
                      .data$p_res,
                      .data$p_rel),
            ## refer to list based on index
            .f = ~ buildDGP(
              # varying parameters
              n = ..1,
              t = ..2,
              strength = ..3,
              rate = ..4,
              balance_dir = ..5,
              balance_res = ..6,
              reliable = ..7,
              export = TRUE,
              patterns = FALSE,
              slopes = TRUE
            ), ## for reproducibility
            .options = furrr::furrr_options(seed = TRUE),
            .progress = TRUE
          )
      )

    # ----------------------------------------------------- #
    # CALCULATE CONFUSION                                   #
    # ----------------------------------------------------- #

    data <- data |>
      tidyr::unnest_wider(col = data) |>

      # step 1: classify
      dplyr::mutate(
        correspondence = purrr::map2(
          .x = .data$data,
          .y = .data$slopes,
          .f =
            ~ dplyr::left_join(.x |>
                                 dplyr::distinct(.data$pid,
                                                 .data$changer),
                               .y,
                               by = "pid") |>
            dplyr::mutate(
              predicted =
                dplyr::case_when(
                  (.data$estimate > 0) &
                    (abs(.data$estimate - 0) >= abs(.data$estimate - strength)) ~
                    1,
                  (.data$estimate < 0)
                  &
                    (abs(.data$estimate - 0) >= abs(.data$estimate + strength)) ~
                    1,
                  .default = 0
                )
            ) |>
            dplyr::select("pid", "changer", "predicted")
        )
      ) |>

      # step 2: confuse
      dplyr::mutate(confusion = purrr::map(
        .x = .data$correspondence,
        .f = ~  tibble::tibble(changer =
                                 c(0, 0, 1, 1),
                               predicted =
                                 c(0, 1, 0, 1)) |>
          dplyr::left_join(
            . |>
              dplyr::count(changer, predicted),
            by = c("changer", "predicted")
          ) |>
          dplyr::mutate(n = ifelse(is.na(.data$n), 0, .data$n))
      )) |>

      # step 3: organize
      dplyr::mutate(confused =
                      purrr::map(.x = .data$confusion,
                                 .f = ~ confusionMatrix(.))) |>
      dplyr::select("sims", "confused") |>
      tidyr::unnest_wider(.data$confused)

    return(data)

  }

# ---------------------------------------------------------------------------- #
