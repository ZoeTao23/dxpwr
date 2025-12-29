#' Sample Size Calculator for Diagnostic Tests
#'
#' Adjust sample size for clustered data in diagnostic accuracy studies.
#'
#' @param n1_unadj (integer): Number of required patients with condition if each patient had only one observation.
#' @param avg_obs_per_patient (numeric): Average number of observations per patient.
#' @param intra_corr (numeric): Average correlation between test results of observations from the same patient.
#' @param R (numeric): Ratio of patients with and without the condition (default: 1).
#' @return Object of class `adjust_by_cluster` containing:
#' \itemize{
#'   \item method (character) - Adjustment method.
#'   \item parameter (list) - Input parameters containing:
#'      \itemize{
#'          \item correlation (numeric): Average correlation between test results of observations from the same patient.
#'          \item avg_observations (numeric): Average number of observations per patient.
#'          \item R (numeric): Ratio of patients with and without the condition.
#'     }
#'   \item sample (list) - Output results containing:
#'      \itemize{
#'          \item n_with_condition (integer) - Unadjusted number of patients with the condition.
#'          \item n_without_condition (integer) - Unadjusted number of patients without the condition.
#'          \item n_unadjusted (integer) - Total unadjusted sample size.
#'          \item n_with_condition (integer) - The anticipated sample size for patients with condition in studies with clustered data.
#'          \item n_without_condition (integer) - The anticipated sample size for patients without condition in studies with clustered data.
#'          \item n_adjusted (integer) - The anticipated total sample size in studies with clustered data.
#'     }
#'   \item call (character) - The original function call.
#' }
#' @importFrom stats qnorm
#' @examples
#' n_adj <- adjust_for_cluster_data(n1_unadj=46, avg_obs_per_patient=1.5, intra_corr=0.5, R=1)
#' print(n_adj)
#' @export
adjust_for_cluster_data <- function(n1_unadj, avg_obs_per_patient, intra_corr, R=1) {


  #  (1) ----- Validate inputs
  if (!is.numeric(n1_unadj) || length(n1_unadj) != 1 || n1_unadj <= 0 || !(n1_unadj == round(n1_unadj))) {
    stop("'n1_unadj' must be an single positive integer value.")
  }

  if (!is.numeric(avg_obs_per_patient) || length(avg_obs_per_patient) != 1 ||  avg_obs_per_patient <= 1) {
    stop("'avg_obs_per_patient' must be a single numeric value in (1, +∞).")
  }

  if (!is.numeric(intra_corr) || length(intra_corr) != 1 || intra_corr <= 0 || intra_corr >= 1) {
    stop("'intra_corr' must be a single numeric value in (0, 1)")
  }

  if (!is.numeric(R) || length(R) != 1 || R <= 0) {
    stop("'R' must be a single positive numeric value.")
  }

  #  (2) ----- Calculation modules
  n0_unadj = ceiling(n1_unadj/R)
  n_unadjusted = n0_unadj + n1_unadj

  obs_required <- n1_unadj * (1 + (avg_obs_per_patient - 1) * intra_corr)

  n1_adj <- ceiling(obs_required / avg_obs_per_patient)
  n0_adj <- ceiling(n1_adj/R)
  n_adjusted <- n1_adj + n0_adj


  #  (3) ----- Structured outputs

  ## method reference
  method <- "Sample size adjusted method for cluster data by Obuchowski (1997)"

  result <-
    structure(
      list(
        method = method,
        parameter = list(
          correlation = intra_corr,
          avg_observations = avg_obs_per_patient,
          R = R
        ),
        sample = list(
          n_unadj_with_condition = n1_unadj,
          n_unadj_without_condition = n0_unadj,
          n_unadjusted = n_unadjusted,
          n_adj_with_condition = n1_adj,
          n_adj_without_condition = n0_adj,
          n_adjusted = n_adjusted
        ),
        call = match.call()
      ),
      class = "adjust_by_cluster"
    )

  return(result)
}
