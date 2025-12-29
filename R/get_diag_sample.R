#' Sample Size Calculator for Diagnostic Tests
#'
#' Calculates the required sample size for:
#'    1. estimating the accuracy of one diagnostic test;
#'    2. compare the accuracies of two diagnostic tests;
#'    3. assess non-inferiority or equivalence of two tests.
#'
#' @param var_function (numeric/vector). Variance function value(s) depending on test type:
#'   - For "one_diagnostic": a single numeric value.
#'   - For "two_diagnostic": a numeric vector of length 2 (V_null and V_alternative).
#'   - For "non_inferiority": a single numeric value.
#'   - For "equivalence": a single numeric value.
#' @param alpha (numeric): Significance level (default: 0.05).
#' @param beta (numeric): Type II error rate.
#' @param delta (numeric): Desired half-width of CI/Minimum detectable difference.
#' @param test_type (character). Type of test:
#'   - `"one_diagnostic"`: To evaluate the diagnostic accuracy of a single test.
#'   - `"two_diagnostic"`: To compare the diagnostic accuracy of two tests.
#'   - `"non_inferiority"`: To assess the non-inferiority of test 1 to test 2.
#'   - `"equivalence"`: To test the equivalence of test 1 and test 2.
#' @return N (integer): Estimated sample size.
#' @importFrom stats qnorm
#' @examples
#' N1 <- get_diag_sample(var_function=0.8, alpha=0.05, beta=0.2, delta=0.1, test_type="one_diagnostic")
#' print(N1)
#' N2 <- get_diag_sample(var_function=c(0.2, 0.3), alpha=0.05, beta=0.2, delta=0.1, test_type="two_diagnostic")
#' print(N2)
#' @references
#' Zhou XH, Obuchowski NA, McClish DK. Statistical Methods in Diagnostic Medicine (2nd ed). Wiley; 2011.
#' @export
get_diag_sample <- function(var_function, alpha=0.05, beta=NULL, delta, test_type){

  # (1) ----- Validate inputs
  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(delta) || length(delta) != 1 || delta <= 0) {
    stop("'delta' must be a single positive numeric value.")
  }

  if (!is.character(test_type) || !(test_type %in% c("one_diagnostic", "two_diagnostic", "equivalence", "inferiority"))) {
    stop("'test_type' must be one of 'one_diagnostic', 'two_diagnostic', 'equivalence', or 'inferiority'.")
  }

  if (test_type != "one_diagnostic" && (is.null(beta) || !is.numeric(beta) || length(beta) != 1|| beta <= 0 || beta >= 1)) {
    stop("'beta' must be a numeric value in (0, 1) for ", test_type, " test.")
  }


  # (2) ----- Calculation module

  ## calculate sample size based on test type
  N <- switch(
    test_type,
    "one_diagnostic" = {
      if (!is.numeric(var_function) || length(var_function) != 1) {
        stop("For 'one_diagnostic', var_function must be a single numeric value.")
      }

      z_alpha <- qnorm(alpha / 2)
      z_beta <- ifelse(is.null(beta), 0, qnorm(beta))
      numerator <- ((z_alpha + z_beta) * sqrt(var_function))^2
      ceiling(numerator / delta^2)
    },

    "two_diagnostic" = {
      if (!is.numeric(var_function) || length(var_function) != 2) {
        stop("For 'two_diagnostic', var_function must be a numeric vector of length 2.")
      }

      z_alpha <- qnorm(alpha / 2)
      z_beta <- qnorm(beta)
      numerator <- (z_alpha * sqrt(var_function[1]) + z_beta * sqrt(var_function[2]))^2
      ceiling(numerator / delta^2)
    },

    "non_inferiority" = {
      if (!is.numeric(var_function) || length(var_function) != 1) {
        stop("For 'non_inferiority', var_function must be a single numeric value.")
      }

      z_alpha <- qnorm(alpha)
      z_beta <- qnorm(beta)
      numerator <- (z_alpha + z_beta)^2 * var_function
      ceiling(numerator / delta^2)
    },

    "equivalence" = {
      if (!is.numeric(var_function) || length(var_function) != 1) {
        stop("For 'equivalence', var_function must be a single numeric value.")
      }

      z_alpha <- qnorm(alpha)
      z_beta <- qnorm(beta)
      numerator <- (z_alpha + z_beta)^2 * var_function
      ceiling(numerator / delta^2)
    },

    stop("Invalid test type. Use 'one_diagnostic', 'two_diagnostic', 'inferiority' or 'equivalence'.")
  )


  #  (3) ----- Structured outputs
  return(N)

}
