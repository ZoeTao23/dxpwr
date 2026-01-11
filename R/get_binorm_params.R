#' Sample Size Calculator for Diagnostic Tests
#'
#' Estimate binormal model parameters (a, b) for Receiver Operating Characteristic (ROC) curve fitting.
#'
#' @param AUC (numeric): Areas under the ROC curve.
#' @param SD_Ratio (numeric): Ratio of standard deviations of the distributions of test results (without vs with condition).
#' @return A list of class `binorm_params` containing:
#' \itemize{
#'   \item a (numeric) - Standardized difference in means of the distributions of test results for patients with and without the condition, (μ_with - μ_without)/σ_with.
#'   \item b (numeric) - Ratio of the standard deviations of the distributions of test results for patients without versus with the condition, σ_without/σ_with.
#'   \item interpretation (list) - The interpretation of the ROC parameters (a, b).
#'   \item call (character) - The original function call.
#' }
#' @importFrom stats qnorm
#' @examples
#' params <- get_binorm_params(AUC=0.8, SD_Ratio=1.2)
#' print(params)
#' @export
get_binorm_params <- function(AUC, SD_Ratio){

  #  (1) ----- validate inputs
  if (!is.numeric(AUC) || length(AUC) != 1 || AUC < 0.5 || AUC > 1) {
    stop("'AUC' must be a single numeric value in [0.5, 1].")
  }

  if (!is.numeric(SD_Ratio) || length(SD_Ratio) != 1 || SD_Ratio <= 0) {
    stop("'SD_Ratio' must be a single positive numeric value.")
  }


  #  (2) ----- calculation modules
  b <- SD_Ratio
  a <- qnorm(AUC) * sqrt(1 + b^2)


  #  (3) ----- Structured outputs
  params <-
    structure(
      list(
        a = a,
        b = b,
        interpretation = list(
          a = "Standardized mean difference: (μ_with - μ_without)/σ_with",
          b = "SD ratio: σ_without/σ_with"
        ),
        call = match.call()
      ),
      class = "binorm_params"
    )

  return(params)

}
