#' Sample Size Calculator for Diagnostic Tests
#'
#' Adjust sample size for prospective diagnostic accuracy studies, accounting for the prevalence of the condition.
#'
#' @param n_unadjusted (integer): The initial unadjusted sample size calculated without considering prevalence.
#' @param prev (numeric): The expected prevalence of the target condition in the population.
#' @param beta (numeric): Type II error rate.
#' @return Object of class `adjust_by_prev` containing:
#' \itemize{
#'   \item method (character) - Adjustment method.
#'   \item parameter (list) - Input parameters containing:
#'      \itemize{
#'          \item power (numeric) - Desired statistical power (1 - Type II error rate).
#'          \item prevalence (numeric) - The expected prevalence of the target condition in the population.
#'     }
#'   \item sample (list) - Output results containing:
#'      \itemize{
#'          \item n_unadjusted (integer) - The initial unadjusted sample size calculated without considering prevalence.
#'          \item n_adjusted (integer) - The total sample size should be recruited.
#'     }
#'   \item call (character) - The original function call.
#' }
#' @importFrom stats qnorm
#' @examples
#' n_adj <- adjust_by_disease_prevalence(sample=350, prev=0.8, beta=0.9)
#' print(n_adj)
#' @export
adjust_by_disease_prevalence <- function(sample, prev, beta) {

  #  (1) ----- Validate inputs
  if (!is.numeric(n_unadjusted) || length(n_unadjusted) != 1 || n_unadjusted <= 0 || !(n_unadjusted == round(n_unadjusted))) {
    stop("'n_unadjusted' must be a single positive integer.")
  }

  if (!is.numeric(prev) || length(prev) != 1 || prev <= 0 || prev >= 1) {
    stop("'prev' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(beta) || length(beta) != 1 || beta <= 0 || beta >= 1) {
    stop("'beta' must be a single numeric value in (0, 1).")
  }


  #  (2) ----- Calculation modules


  ## calculate critical value
  z_beta <- qnorm(beta)

  ## calculate quadratic equation coefficients
  a <- prev^2
  b <- -2 * prev * sample - prev * (1 - prev) * z_beta^2
  c <- sample^2

  ## adjust sample size
  n_adjusted <- ceiling((-b + sqrt(b^2 - 4 * a * c)) / (2 * a))


  #  (3) ----- Structured outputs

  ## method reference
  method <- "Sample size adjusted method based on prevalence"

  result <-
    structure(
      list(
        method = method,
        parameter = list(
          power = 1 - beta,
          prevalence = prev
        ),
        sample = list(
          n_unadjusted = n_unadjusted,
          n_adjusted = n_adjusted
        ),
        call = match.call()
      ),
      class = "adjust_by_prev"
    )

  return(result)

}
