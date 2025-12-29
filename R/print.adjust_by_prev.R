#' Print method
#'
#' Display a human-readable summary of adjust_by_prev objects containing sample size before and after adjustment.
#'
#' @param x A `adjust_by_prev` object returned by [adjust_by_disease_prevalence()].
#' @return Invisibly returns for `adjust_by_prev` object.
#' @examples
#' n_adj <- adjust_by_disease_prevalence(n_unadjusted=350, prev=0.8, beta=0.9)
#' print(n_adj)
#' n_adj  # same as print(n_adj)
#' @export
print.adjust_by_prev <- function(x) {

  cat("\033[1mAdjustment Method:\033[0m ", x$method, "\n\n")

  cat("\033[1mCall:\033[0m\n")
  print(x$call)
  cat("\n")

  cat("\033[1mParameters:\033[0m\n")
  cat(sprintf("  %-10s: %.*f\n", "Power",2, x$parameter$power))
  cat(sprintf("  %-10s: %.*f\n", "Prevalence", 2, x$parameter$prevalence))

  cat("\n\033[1mSample Size:\033[0m\n")
  cat(sprintf("  %-10s: n = %d\n", "Before", x$sample$n_unadjusted))
  cat(sprintf("  %-10s: n = %d\n", "After", x$sample$n_adjusted))

  invisible(x)
}
