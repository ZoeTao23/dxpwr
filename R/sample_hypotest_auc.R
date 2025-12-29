#' Sample Size Calculator for Diagnostic Tests
#'
#'Calculate the required sample size for testing the hypothesis that the Area Under the ROC Curve is equal to a particular value.
#'
#' @title Sample Size Calculation for Fixed AUC Comparison
#' @description Computes the required sample size for testing whether the AUC is equal to a specific value.
#'
#' @param A_null Numeric, assumed null hypothesis AUC (0.5 < A_null < 1).
#' @param A_alter Numeric, alternative hypothesis AUC (0.5 < A_alter < 1).
#' @param b Numeric, standard deviation ratio (default is 1, must be > 0).
#' @param alpha Numeric, significance level (0 < alpha < 1).
#' @param beta Numeric, Type II error (0 < beta < 1, power = 1 - beta).
#' @param alternative Character, type of test, must be one of "less", "greater", or "two.sided" (default: "two.sided").
#' @param R Numeric, ratio of non-diseased to diseased subjects (default is 1, must be > 0).
#' @return An object of class "diag_sample_size_html" containing:
#'   - `sample_size` (list): Required sample size.
#'   - `parameters` (list): Input parameters.
#'   - `html_report` (html): Sample size calculation report.
#' @import stats
#' @examples
#' n <- sample_hypotest_auc(A_null=0.5, A_alter=0.8, alpha=0.05, beta=0.2, alternative="two.sided", R=1)
#' print(n)
#' @export
sample_hypotest_auc <- function(A_null, A_alter, b=1, alpha=0.05, beta, alternative="two.sided", R=1) {

  #--------------------------#
  #  Sample size calculation when testing hypothesis that auc is equal to particular value

  #  Args:
  #     A_null (numeric): AUC under the null hypothesis
  #     A_alter (numeric): AUC under the alternative hypothesis
  #     b (numeric): the ratio of standard deviations (Sx/Sy)
  #     alpha (numeric): Significance level (default: 0.05)
  #     beta (numeric): Type II error rate
  #     alternative (character): alternative hypothesis, must be one of "less", "greater" or"two-sided"
  #     R (numeric): the ratio of patients with and without the condition

  #  Returns:
  #     sample size: (with condition, without condition, total)

  #--------------------------#

  #  (1) ----- validate inputs

  if (!is.numeric(A_null) || length(A_null) != 1 || A_null <= 0 || A_null >= 1) {
    stop("'A_null' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(A_alter) || length(A_alter) != 1  || A_alter <= 0 || A_alter >= 1) {
    stop("'A_alter' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(b) || length(b) != 1  || b <= 0) {
    stop("'b' must be a single positive numeric value.")
  }

  if (!is.numeric(alpha) || length(alpha) != 1  || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(beta) || length(beta) != 1  || beta <= 0 || beta >= 1) {
    stop("'beta' must be a single numeric value in (0, 1).")
  }

  if (!is.character(alternative) || length(alternative) != 1  || !(alternative %in% c("two.sided", "greater", "less"))) {
    stop("'alternative' must be one of 'two.sided', 'greater' or 'less'.")
  }

  if (!is.numeric(R) || length(R) != 1  || R <= 0) {
    stop("'R' must be a single positive numeric value.")
  }

  #  (2) ----- calculation modules

  ## calculate variance function
  a0 <- get_binorm_params(A_null, b)$a
  V0_theta <- (0.0099) * exp(-a0^2 / 2)*((5 * a0^2 + 8)+(a0^2 + 8) / R)

  aA <- get_binorm_params(A_alter, b)$a
  VA_theta <- (0.0099) * exp(-aA^2 / 2)*((5 * aA^2 + 8)+(aA^2 + 8) / R)

  var_function <- c(V0_theta, VA_theta)
  delta <- abs(A_null - A_alter)

  ## adjust alpha for alternative test
  alpha <- switch(
    alternative,
    "less" = 2 * alpha,
    "greater" = 2 * alpha,
    "two.sided" = alpha,
    stop("Invalid 'alternative' value. Use 'less', 'greater', or 'two.sided'.")
  )

  ## calculate sample size
  n <- get_diag_sample(
    var_function = var_function,
    alpha = alpha,
    beta = beta,
    delta = delta,
    test_type = "multi_diagnostic"
  )

  ## calculate patients with and without condition
  n_with_condition <- ceiling(n)
  n_without_condition <- ceiling(n*R)
  N_total <- n_with_condition + n_without_condition


  #  (3) ----- structured outputs

  ## method
  method <- switch(
    alternative,
    "less" = paste0("Non-inferiority test: AUC < ",A_null),
    "greater" = paste0("Superiority test: AUC > ",A_null),
    "two.sided" = paste0("Equivalence test: AUC = ",A_null)
  )

  structure(
    list(
      alpha = alpha,
      power = ifelse(is.null(beta),"NULL", 1 - beta),
      method = method,
      sample_size = N_total,
      n_with_condition = n_with_condition,
      n_without_condition = n_without_condition
    ),
    class = "Single diagnostic test")

}
