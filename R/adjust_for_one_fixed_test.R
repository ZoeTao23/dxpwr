#' Sample Size Calculator for Diagnostic Tests
#'
#' Adjust sample size for an unpaired study design comparing two diagnostic tests, when sample size for one test is fixed.
#'
#' @param n_known_test (integer): The known/predetermined sample size for one test, default as test 1.
#' @param n_unadjusted (integer): The desired total sample size for the study.
#' @return Object of class `adjust_by_fixed_size` containing:
#' \itemize{
#'   \item method (character) - Adjustment method.
#'   \item parameter (list) - Input parameters containing:
#'      \itemize{
#'          \item n_known_test (integer) - The known/predetermined sample size for one test.
#'          \item n_unadjusted (integer) - The desired total sample size for the study.
#'     }
#'   \item sample (list) - Output results containing:
#'      \itemize{
#'          \item n_unadj_known_test (integer) - The known/predetermined sample size for one test.
#'          \item n_unadj_unknown_test (integer) - Unadjusted sample size for another test.
#'          \item n_unadjusted (integer) - The desired total sample size for the study.
#'          \item n_adj_known_test (integer) - The known/predetermined sample size for one test.
#'          \item n_adj_unknown_test (integer) - The desired sample size for another test.
#'          \item n_adjusted (integer) - The desired total sample size after adjustment.
#'     }
#'   \item call (character) - The original function call.
#' }
#' @importFrom stats qnorm
#' @examples
#' n_adj <- adjust_for_one_fixed_test(n_known_test=52, n_unadjusted=77)
#' print(n_adj)
#' @export
adjust_for_one_fixed_test <- function(n_known_test, n_unadjusted) {

  #  (1) ----- Validate inputs
  if (!is.numeric(n_known_test) || length(n_known_test) != 1 || n_known_test <= 0 || !(n_known_test == round(n_known_test))) {
    stop("'n_known_test' must be an single positive integer value.")
  }

  if (!is.numeric(n_unadjusted) || length(n_unadjusted) != 1 || n_unadjusted <= 0 || !(n_unadjusted == round(n_unadjusted)) || n_unadjusted <= n_known_test) {
    stop("'n_unadjusted' must be an single integer value larger than 'n_known_test'.")
  }


  #  (2) ----- Calculation modules
  n_unknown_test <- ceiling((n_unadjusted*n_known_test) / (2*n_known_test - n_unadjusted) )
  n_adjusted <- n_known_test + n_unknown_test


  #  (3) ----- Structured outputs

  ## reference method
  method <- "Sample size adjusted method in unpaired design when sample size for one test is fixed by Cohen (1977)"

  structure(
    list(
      method = method,
      parameter = list(
        n_known_test = n_known_test,
        n_unadjusted = n_unadjusted
      ),
      sample = list(
        n_unadj_known_test = n_known_test,
        n_unadj_unknown_test = n_unadjusted - n_known_test,
        n_unadjusted = n_unadjusted,
        n_adj_known_test = n_known_test,
        n_adj_unknown_test = n_unknown_test,
        n_adjusted = n_adjusted
      ),
      call = match.call()
    ),
    class = "adjust_by_fixed_size"
  )

}



