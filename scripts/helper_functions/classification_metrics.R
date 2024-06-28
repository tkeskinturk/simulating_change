
# - promises and pitfalls of panel change ------------------------------------ #
# - helper: classification metrics ------------------------------------------- #

### note: this script
###       (a) takes a confusion matrix of true/simulated counts as its input,
###       (b) calculates kappa, sensitivity, and specificity.

require(tidyverse)

# FUNCTION ------------------------------------------------------------------- #

confused_outcomes <-
  
  function(confusion_matrix) {
    
    # --- calculate kappa
    
    # accuracy
    p_0 =
      (confusion_matrix[4, 3] 
       + 
       confusion_matrix[1, 3]
      ) / 
      sum(confusion_matrix$n)
    
    # expected classification
    p_e =
      (
        (
        (confusion_matrix[1, 3] + confusion_matrix[2, 3])
        *
        (confusion_matrix[1, 3] + confusion_matrix[3, 3])
        )
      +
        (
        (confusion_matrix[3, 3] + confusion_matrix[4, 3])
        *
        (confusion_matrix[2, 3] + confusion_matrix[4, 3])
        )
      ) /
      (sum(confusion_matrix$n) ^ 2)
    
    # kappa score
    kappa = (p_0 - p_e) / (1 - p_e)
    
    # --- calculate sensitivity
    sensitivity =
      confusion_matrix[4, 3] /
      (confusion_matrix[4, 3] + confusion_matrix[3, 3])
    
    # --- calculate specificity
    specificity =
      confusion_matrix[1, 3] /
      (confusion_matrix[1, 3] + confusion_matrix[2, 3])
    
    # --- export
    output <- list(kappa = kappa$n,
                   sensitivity = sensitivity$n,
                   specificity = specificity$n)
    
    return(output)
  }

# ---------------------------------------------------------------------------- #
