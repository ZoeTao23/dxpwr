#' Sample Size Calculator for Diagnostic Tests
#'
#' Estimate required sample size to identify a suitable cutoff value in diagnostic accuracy studies.
#'
#' @param spec_min (numeric): The acceptable minimum specificity in order for the cutoff to be useful.
#' @param sens_at_spec_min (numeric): The conjectured sensitivity at the acceptable minimum specificity.
#' @param sens_min (numeric): The acceptable minimum sensitivity in order for the cutoff to be useful.
#' @param alpha (numeric): Significance level (default: 0.05).
#' @param beta (numeric): Type II error rate.
#' @param b (numeric): Ratio of the standard deviations of the distributions of test results for patients without versus with the condition (default: 1).
#' @return Object of class `optimize_for_cutoff` containing:
#' \itemize{
#'   \item method (character) - Adjustment method.
#'   \item parameter (list) - Input parameters containing:
#'      \itemize{
#'          \item min_spec (numeric) - The acceptable minimum specificity in order for the cutoff to be useful.
#'          \item sens_at_min_spec (numeric) - The conjectured sensitivity at the acceptable minimum specificity.
#'          \item sens_min (numeric) - The acceptable minimum sensitivity in order for the cutoff to be useful.
#'          \item alpha (numeric) - Significance level.
#'          \item power (numeric) - Desired statistical power (1 - Type II error rate).
#'          \item b (numeric) - Ratio of the standard deviations of the distributions of test results for patients without versus with the condition, σ_without/σ_with.
#'     }
#'   \item sample (list) - Output results containing:
#'      \itemize{
#'          \item n_total (integer) - The required total sample size.
#'     }
#'   \item call (character) - The original function call.
#' }
#' @importFrom stats qnorm
#' @examples
#' n_total <- optimize_for_threshold(spec_min=0.9, sens_at_spec_min=0.6, sens_min=0.5, alpha=0.05, beta=0.2, b=1)
#' print(n_total)
#' @export
optimize_for_threshold <- function(spec_min, sens_at_spec_min, sens_min, alpha, beta, b=1) {

  #  (1) ----- Validate inputs
  if (!is.numeric(spec_min) || length(spec_min) != 1 || spec_min <= 0 || spec_min >= 1) {
    stop("'spec_min' must be a single numeric value in (0, 1)")
  }

  if (!is.numeric(sens_at_spec_min) || length(sens_at_spec_min) != 1 || sens_at_spec_min <= 0 || sens_at_spec_min >= 1) {
    stop("'sens_at_spec_min' must be a single numeric value in (0, 1)")
  }

  if (!is.numeric(sens_min) || length(sens_min) != 1 || sens_min <= 0 || sens_min >= 1) {
    stop("'sens_min' must be a single numeric value in (0, 1)")
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value in (0, 1)")
  }

  if (!is.numeric(beta) || length(beta) != 1 || beta <= 0 || beta >= 1) {
    stop("'beta' must be a single numeric value in (0, 1)")
  }

  if (!is.numeric(b) || length(b) != 1 || b <= 0) {
    stop("'b' must be a single positive numeric value.")
  }


  #  (2) ----- Calculation modules

  ## calculate critical value
  z_alpha <- qnorm(sqrt(1 - alpha))
  z_beta <- qnorm(1 - beta)

  ## calculate sample size
  lambda <- qnorm(sens_at_spec_min) - qnorm(sens_min)
  v_x <- b * sqrt(1 + 1 / 2 * qnorm(spec_min)^2)
  v_y <- sqrt(1 + 1 / 2 *qnorm(sens_min)^2)

  n_total <- ((sqrt(2) * z_alpha + z_beta)^2 * (v_x + v_y)^2) / lambda^2
  n_total <- ceiling(n_total)


  #  (3) ----- Structured outputs

  ## reference method
  method <- "Sample size method for determining a suitable cutoff value by Schäfer (1989)"

  result <- structure(
    list(
      method = method,
      parameter = list(
        min_spec = spec_min,
        sens_at_min_spec = sens_at_spec_min,
        min_sens = sens_min,
        alpha = alpha,
        power = 1 - beta,
        b = b
      ),
      sample = list(
        n_total = n_total
      ),
      call = match.call()
    ),
    class = "optimize_for_threshold"
  )

  return(result)
}



