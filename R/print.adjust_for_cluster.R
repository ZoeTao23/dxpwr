#' Print method
#'
#' Display a human-readable summary of adjust_by_cluster objects containing sample size before and after adjustment.
#'
#' @param x A `adjust_by_cluster` object returned by [adjust_for_cluster_data()].
#' @return Invisibly returns for `adjust_by_cluster` object.
#' @examples
#' n_adj <- adjust_for_cluster_data(n1_unadj=46, avg_obs_per_patient=1.5, intra_corr=0.5)
#' print(n_adj)
#' n_adj  # same as print(n_adj)
#' @export
print.adjust_by_cluster <- function(x) {

  cat("\033[1mAdjustment Method:\033[0m ", x$method, "\n\n")

  cat("\033[1mCall:\033[0m\n")
  print(x$call)
  cat("\n")

  cat("\033[1mParameters:\033[0m\n")
  cat(sprintf("  %-22s: %.*f\n", "Correlation", 2, x$parameter$correlation))
  cat(sprintf("  %-22s: %.*f\n", "Average observations",2, x$parameter$avg_observations))
  cat(sprintf("  %-22s: %.*f\n", "R", 2, x$parameter$R))

  cat("\n\033[1mSample Size:\033[0m\n")
  cat(sprintf("  %-6s: n_with = %d, n_without = %d, total = %d\n", "Before", x$sample$n_unadj_with_condition, x$sample$n_unadj_without_condition, x$sample$n_unadjusted))
  cat(sprintf("  %-6s: n_with = %d, n_without = %d, total = %d\n", "After", x$sample$n_adj_with_condition, x$sample$n_adj_without_condition, x$sample$n_adjusted))

  invisible(x)
}
