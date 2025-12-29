#' Sample Size Calculator for Diagnostic Tests
#'
#' Calculate the required sample size for comparing the positive predictive value (PPV) or negative predictive value (NPV) between two diagnostic tests.
#'
#' @param theta2 Numeric. NPV or PPV of test 2 (0 < theta2 < 1).
#' @param p Numeric vector (length = 8). Probability components required for calculation.
#' \describe{
#'     For patients without the condition：
#'     \item{p[1]}{the proportion of patients testing positive on both tests}
#'     \item{p[2]}{the proportion testing positive on test 1 and negative on test 2}
#'     \item{p[3]}{the proportion testing positive on test 2 and negative on test 1}
#'     \item{p[4]}{the proportion testing negative on both tests}
#'     For patients with the condition：
#'     \item{p[5]}{the proportion of patients testing positive on both tests}
#'     \item{p[6]}{the proportion testing positive on test 1 and negative on test 2}
#'     \item{p[7]}{the proportion testing positive on test 2 and negative on test 1}
#'     \item{p[8]}{the proportion testing negative on both tests}
#'   }
#' @param alpha Numeric. Significance level (0 < alpha < 1).
#' @param beta Numeric. Type II error rate (0 < beta < 1).
#' @param gamma Numeric. A specific value of interest for rPPV under the alternative hypothesis (gamma > 0).
#' @param delta Numeric. The null hypothesis is rPPV = delta. (In most cases, let delta = 1).
#' @param paired Logical. Whether the design is paired (`TRUE`) or unpaired (`FALSE`).
#' @param choice Character. `"PPV"` for positive predictive value, `"NPV"` for negative predictive value.
#' @param alternative Character. `"less"`, `"greater"`, or `"two.sided"`.
#' @return An object of class "diag_sample_size_html" containing:
#'   - `sample_size` (list): Required sample size.
#'   - `parameters` (list): Input parameters.
#'   - `html_report` (html): Sample size calculation report.
#'
#' @examples
#' n1 <- sample_compare_ppvnpv(theta2=0.8, p=c(0.4,0,0.52,0.44,0.7,0,0.1,0.2), alpha=0.05, beta=0.2, delta=1, gamma=1.1, paired=TRUE, metric="PPV", alternative="two.sided")
#' print(n1)
#' n2 <- sample_compare_ppvnpv(theta2=0.8, p=c(0.4,0,0.52,0.44,0.7,0,0.1,0.2), alpha=0.10, beta=0.2, delta=1, gamma=1.37, paired=TRUE, metric="NPV", alternative="two.sided")
#' print(n2)
sample_compare_ppvnpv <- function(metric="PPV", theta2, p, alpha=0.05, beta, gamma, delta, paired=TRUE, alternative="two.sided") {

  #--------------------------#
  #  Sample size calculation when comparing test's positive/negative predictive values

  #  Args:
  #     metric (character): index to be compared between two diagnostic tests, must be either "PPV" or "NPV"
  #     theta2 (numeric): the expected NPV (or PPV) of test2
  #     alpha (numeric): Significance level (default: 0.05)
  #     beta (numeric): Type II error rate
  #     p (vector): the proportion vector with 8 elements:
  #       for patients without the condition,
  #         p1: the pr of patients (test1 +, test2 +);
  #         p2: the pr of patients (test1 +, test2 -);
  #         p3: the pr of patients (test1 -, test2 +);
  #         p4: the pr of patients (test1 -, test2 -).
  #       for patients with the condition,
  #         p5: the pr of patients (test1 +, test2 +);
  #         p6: the pr of patients (test1 +, test2 -);
  #         p7: the pr of patients (test1 -, test2 +);
  #         p8: the pr of patients (test1 -, test2 -).
  #     gamma (numeric): a specific value of interest for rPPV/rNPV under the alternative hypothesis (gamma > 0)
  #     delta (numeric): a specific value of interest in the null hypothesis (In most cases, let delta = 1).
  #     paired (logical): paired study design or unpaired study design, must be "True" or "False"


  #  Returns:
  #     sample size: (test1, test2, total)

  #--------------------------#

  #  (1) ----- validate inputs

  if (length(metric) != 1 || !metric %in% c("PPV", "NPV")) {
    stop("'metric' must be either 'PPV' or 'NPV'.")
  }

  if (!is.numeric(theta2) || length(theta2) != 1 || theta2 <= 0 || theta2 >= 1) {
    stop("'theta2' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(p) || length(p) != 8 || any(p < 0 | p > 1)) {
    stop("'p' must be a numeric vector of length 8 with values in [0, 1].")
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(beta) || length(beta) != 1 || beta <= 0 || beta >= 1) {
    stop("'beta' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(gamma) || length(gamma) != 1 || gamma <= 0) {
    stop("'gamma' must be a single positive numeric value.")
  }

  if (!is.numeric(delta) || length(delta) != 1 || delta <= 0) {
    stop("'delta' must be a single positive numeric value.")
  }

  if (!is.logical(paired) || length(paired) != 1) {
    stop("'paired' must be a single logical value TRUE or FALSE.")
  }

  if (length(alternative) != 1 || !alternative %in% c("less", "greater", "two.sided")) {
    stop("'alternative' must be one of 'less', 'greater', or 'two.sided'.")
  }


  #  (2) -----  calculation modules
  null_hypothesis <- switch(
    alternative,
    "less" = paste0("null hypothesis: r", metric," > ", delta),
    "greater" = paste0("null hypothesis: r", metric, " < ", delta),
    "two.sided" = paste0("null hypothesis: r", metric, " = ", delta)
  )
  message(print(null_hypothesis))


  ## adjust alpha for one-sided or two-sided test
  alpha <- switch(
    alternative,
    "less" = 2 * alpha,
    "greater" = 2 * alpha,
    "two.sided" = alpha,
    stop("Invalid 'alternative' value. Use 'less', 'greater', or 'two.sided'.")
  )

  ## calculate critical values
  z_alpha <- qnorm(1 - alpha / 2)
  z_beta <- qnorm(1 - beta)

  ## calculate sample size
  if(metric == "PPV") {

    m1 <- ((z_alpha + z_beta) / log10(gamma / delta))^2
    m2 <- ((p[5] + p[6]) * (p[5] + p[7])) ^ (-1)
    m3 <- 2 * (p[7] + p[3]) * gamma * theta2^2 +
      (p[5] * (1 - gamma) - p[6]) * theta2 +
      p[6] + p[7] * (1 - 3 * gamma * theta2)

  } else if(metric == "NPV") {

    m1 <- ((z_alpha + z_beta)^2 / log10(gamma / delta))^2
    m2 <- ((p[2] + p[4]) * (p[3] + p[4])) ^ (-1)
    m3 <- -2 * (p[4] + p[8]) * gamma * theta2^2 +
      (-p[3] + p[4] - gamma * (p[2] - p[4])) * theta2 +
      p[2] + p[3]
  } else {
    stop("Invalid 'choice' value. Use 'PPV' or 'NPV'.")
  }

  ## calculate patients in each test
  n <- ceiling(m1*m2*m3)
  n_test1 <- n_test2 <- ceiling(n/2)
  N_total <- n_test1 + n_test2


  #  (3) ----- structured outputs

  ## method reference
  method <- "Sample size method for comparing test's positive/negative predictive values by Moskowize and Pepe(2006)"

  structure(
    list(
      alpha = alpha,
      power = ifelse(is.null(beta),"NULL",1-beta),
      method = paste0(method, ifelse(paired, "- Paired", "- Unpaired")),
      sample_size = N_total,
      n_test1 = n_test1,
      n_test2 = n_test2
    ),
    class = "Two diagnostic tests"
  )

}
