#' Sample Size Calculator for Diagnostic Tests
#'
#' Calculate the required sample size for testing whether the Area Under the ROC Curve (AUC) of a new test (test1) is non-inferior to an existing test (test2).
#'
#' @param A1_alter (numeric): The expected auc of test1
#' @param A2_alter (numeric): The expected auc of test2
#' @param r (numeric): The correlation between the tests because of the paired design
#' @param delta (numeric): The smallest unacceptable difference in accuracy
#' @param alpha (numeric): Significance level (default: 0.05)
#' @param beta (numeric): Type II error rate
#' @param paired (character): Paired study design or unpaired study design, must be "True" or "False"
#' @param R (numeric): Ratio of patients with and without the condition
#' @return An object of class "diag_sample_size_html" containing:
#'   - `sample_size` (list): Required sample size.
#'   - `parameters` (list): Input parameters.
#'   - `html_report` (html): Sample size calculation report.
#' @examples
#' n <- sample_non_inferiority_auc(A1=0.9, A2=0.88, rN=0.5, rD=0.5, delta=0.05, alpha=0.05, beta=0.2, R=1, paired=TRUE, dist="binorm")
#' print(n)
#' @export
sample_non_inferiority_auc <- function(A1,b1=1,A2,b2=1,r=NULL, rN=NULL, rD=NULL, delta, alpha=0.05, beta, R=1, paired=TRUE, dist="any") {

  #--------------------------#
  #  Sample size calculation when testing the auc of a new test (test1) is non-inferior to an existing test (test2)

  #  Args:
  #     A1_alter (numeric): the expected auc of test1
  #     A2_alter (numeric): the expected auc of test1
  #     r (numeric): the correlation between the tests because of the paired design
  #     delta (numeric): the smallest unacceptable difference in accuracy
  #     alpha (numeric): Significance level (default: 0.05)
  #     beta (numeric): Type II error rate
  #     paired (character): paired study design or unpaired study design, must be "True" or "False"
  #     R (numeric): the ratio of patients with and without the condition

  #  Returns:
  #     sample size: (test1 with condition, test1 without condition, test2 with condition, test 2 without condition, total)
  #--------------------------#


  #  (1) ----- validate inputs
  if (!is.numeric(A1) || length(A1) != 1 || A1 <= 0 || A1 >= 1) {
    stop("'A1' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(A2) || length(A2) != 1 || A2 <= 0 || A2 >= 1) {
    stop("'A2' must be a single numeric value in (0, 1).")
  }

  if (!is.logical(paired) || length(paired) != 1) {
    stop("'paired' must be a single logical value TRUE or FALSE.")
  }

  if(dist == "any") {

    if (!is.numeric(r) || length(r) != 1 || r < -1 || r > 1) {
      stop("'r' must be a single numeric value in [-1, 1].")
    }

  } else if(dist == "binorm") {

    if (!is.numeric(b1) || length(b1) != 1 || b1 <= 0) {
      stop("'b1' must be a single positive numeric value.")
    }

    if (!is.numeric(b2) || length(b2) != 1 || b2 <= 0) {
      stop("'b2' must be a single positive numeric value.")
    }

    if (!is.numeric(rN) || length(rN) != 1  || rN < -1 || rN > 1) {
      stop("'rN' must be a single numeric value in [-1, 1].")
    }

    if (!is.numeric(rD) || length(rD) != 1 || rD < -1 || rD > 1) {
      stop("'rD' must be a single numeric value in [-1, 1].")
    }

  } else{

    stop("'dist' must be either 'any' or 'binorm'.")
  }

  if (!is.numeric(delta) || length(delta) != 1 || delta < 0 || delta > 1) {
    stop("'delta' must be a single numeric value in [0, 1].")
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(beta) || length(beta) != 1 || beta <= 0 || beta >= 1) {
    stop("'beta' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(R) || length(R) != 1 || R <= 0) {
    stop("'R' must be a single positive numeric value.")
  }


  #  (2) ----- calculation modules

  message("null_hypothesis: AUC_test1 - AUC_test2 ≥ ", delta)
  message("alternative_hypothesis: AUC_test1 - AUC_test2 < ", delta)

  ## define correlation coefficients for paired and unpaired designs
  if(paired){

    rN <- ifelse(is.null(NULL), 0.5, rN)
    rD <- ifelse(is.null(NULL), 0.5, rD)
    r <- ifelse(is.null(NULL), 0.5, r)

  } else {
    r <- rN <- rD <- 0
  }

  ## calculate variance function
  if(dist == "any") {

    VA_delta_A <- A1 * (1 - A1) + A2 * (1 - A2) - 2 * r * sqrt(A1 * (1-A1) * A2 * (1-A2))
    var_function <- VA_delta_A

  } else {

    a1 <- get_binorm_params(A1, b1)$a
    a2 <- get_binorm_params(A2, b2)$a

    VA <- function(a) {
      V <- (0.0099) * exp(-a^2 / 2) * ((5 * a^2 + 8) + (a^2 + 8) / R)
      return(V)
    }

    VA_A1 <- VA(a1)
    VA_A2 <- VA(a2)


    C_delta_A <- function(a1, a2, rD, rN, R) {

      term1 <- exp(-(a1^2 + a2^2) / 4) / 12.5664 * (rD + rN / R + rD^2* a1 * a2 / 2)
      term2 <- exp(-(a1^2 + a2^2) / 4) / 50.2655 * a1 * a2 * (rN^2 + R* rD^2) / 2 / R
      term3 <- - exp(-(a1^2+a2^2) / 4) / 25.1327 * rD^2 * a1 * a2

      C <- term1 + term2 + term3
      return(C)
    }
    CA_delta_A <- ifelse(paired, C_delta_A(a1, a2, rD, rN, R), 0)

    VA_delta_A <- VA_A1 + VA_A2 - 2 * CA_delta_A

    var_function <- VA_delta_A


  }

  ## calculate delta
  delta <- abs(A1 - A2 - delta)

  ## calculate sample size
  n <- get_diag_sample(
    var_function = var_function,
    alpha = alpha,
    beta = beta,
    delta = delta,
    test_type = "non_inferiority"
  )

  ## calculate sample size with and without condition
  if (paired) {
    n_test1_with_condition <- n_test2_with_condition <- ceiling(n)
    n_test1_without_condition <- n_test2_without_condition <- ceiling(R*n)
    N_total <- n_test1_with_condition + n_test1_without_condition

  } else {
    n_test1_with_condition <- n_test2_with_condition <- ceiling(n)
    n_test1_without_condition <- n_test2_without_condition <- ceiling(R*n)
    N_total <- 2*(n_test1_with_condition + n_test1_without_condition)
  }


  #  (3) ----- structured outputs
  method <- "Sample size method for testing whether auc of a new test is non-inferior to an existing test"

  structure(
    list(
      alpha = alpha,
      power = ifelse(is.null(beta),"NULL",1-beta),
      method = paste0(method, "- Unpaired"),
      sample_size = N_total,
      n_test1_with_condition = n_test1_with_condition,
      n_test1_without_condition = n_test1_without_condition,
      n_test2_with_condition = n_test2_with_condition,
      n_test2_without_condition = n_test2_without_condition
    ),
    class = "Non-inferiority test")


}
