#' Sample Size Calculator for Diagnostic Tests
#'
#' Calculate the required sample size for testing whether the sensitivity/specificity of a new test (test1) is non-inferior to an existing test (test2).
#'
#' @param theta1 (numeric): The expected sensitivity (or specificity) of a new test (test 1).
#' @param theta2 (numeric): The expected sensitivity (or specificity) of an existing test (test 2).
#' @param prob_t1pos_give_t2pos (numeric): Conditional probability of test 1 being positive when test 2 being positive.
#' @param alpha (numeric): Significance level (default: 0.05).
#' @param beta (numeric): Type II error rate.
#' @param delta (numeric): The smallest difference in accuracy which is not acceptable.
#' @param paired (logical): Whether the design is paired (`TRUE`) or unpaired (`FALSE`).
#' @return An object of class "diag_sample_size_html" containing:
#'   - `sample_size` (list): Required sample size.
#'   - `parameters` (list): Input parameters.
#'   - `html_report` (html): Sample size calculation report.
#' @examples
#' n <- sample_non_inferior_sesp(theta1=0.8, theta2=0.9, prob_t1pos_give_t2pos=0, delta=0.1, alpha=0.05, beta=0.2, R=1, paired=TRUE)
#' print(n)
#' @export
sample_non_inferior_sesp <- function(theta1, theta2, prob_t1pos_give_t2pos=NULL, delta, alpha=0.05, beta, R=1, paired=TRUE) {


  #--------------------------#
  #  Sample size calculation when testing the sensitivity/specificity of a new test (test1) is non-inferior to an existing test (test2)

  #  Args:
  #     theta1 (numeric): the expected sensitivity (or specificity) of test1
  #     theta2 (numeric): the expected sensitivity (or specificity) of test2
  #     prob_t1pos_give_t2pos (numeric): probability that test 1 is positive given that test 2 is positive
  #     alpha (numeric): Significance level (default: 0.05)
  #     beta (numeric): Type II error rate
  #     delta (numeric): the smallest unacceptable difference in accuracy
  #     paired (logical): paired study design or unpaired study design, must be "True" or "False"


  #  Returns:
  #     sample size: (test1, test2, total)
  #--------------------------#


  #  (1) ----- validate inputs

  if (!is.numeric(theta1) || length(theta1) != 1 || theta1 <= 0 ||theta1 >= 1) {
    stop("theta1 must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(theta2) || length(theta2) != 1 || theta2 <= 0 || theta2 >= 1) {
    stop("theta2 must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(prob_t1pos_give_t2pos) || length(prob_t1pos_give_t2pos) != 1 ||  prob_t1pos_give_t2pos < 0 || prob_t1pos_give_t2pos > 1) {
    stop("prob_t1pos_give_t2pos must be a single numeric value in [0, 1].")
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("alpha must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(beta) || length(beta) != 1 || beta <= 0 || beta >= 1) {
    stop("beta must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(delta) || length(delta) != 1 || delta < 0 || delta > 1) {
    stop("delta must be a single numeric value in [0, 1].")
  }

  if (!is.logical(paired) || length(paired) != 1) {
    stop("'paired' must be a single logical value TRUE or FALSE.")
  }


  #  (2) -----  calculation modules

  message(paste0("null_hypothesis: theta_test1 - theta_test2 ≥ ", delta))
  message(paste0("alternative_hypothesis: theta_test1 - theta_test2 < ", delta))

  ## calculate variance fucntion
  V0_delta_theta <- theta1 + theta2 - 2 * theta2 * prob_t1pos_give_t2pos
  VA_delta_theta <- V0_delta_theta - (theta1 - theta2)^2

  var_function <- VA_delta_theta

  ## calculate delta
  delta <- abs(theta1 - theta2 - delta)

  ## calculate sample size
  n <- get_diag_sample(
    var_function = var_function,
    alpha = alpha,
    beta = beta,
    delta = delta,
    test_type = "non_inferiority")

  ## calculate number of patients for each test
  if (paired) {
    n_with_condition <- ceiling(n)
    n_without_condition <- ceiling(R*n)
    N_total <- n_with_condition + n_without_condition

  } else {
    n_test1_with_condition <- n_test2_with_condition <- ceiling(n)
    n_test1_without_condition <- n_test2_without_condition <- ceiling(R*n)
    N_total <- 2*(n_test1_with_condition + n_test1_without_condition)
  }


  #  (3) ----- structured outputs

  ## method reference
  method <- "Sample size method for testing whether sensitivity/specificity of a new test is non-inferior to an existing test"

  if (paired) {

    structure(
      list(
        alpha = alpha,
        power = ifelse(is.null(beta), "NULL", 1-beta),
        method = paste0(method, "- Paired"),
        sample_size = N_total,
        n_test1 = n_with_condition,
        n_test2 = n_without_condition
      ),
      class = "Non-inferiority test")

  } else{

    structure(
      list(
        alpha = alpha,
        power = ifelse(is.null(beta),"NULL",1-beta),
        method = paste0(method, "- Unpaired"),
        sample_size = N_total,
        n_test1_with_condition = n_test1_with_condition,
        n_test1_without_condition = n_test1_without_condition,
        n_test2_with_condition = n_test2_with_condition,
        n_test2_without_condition = n_test2_without_condition,
      ),
      class = "Non-inferiority test")
  }

}
