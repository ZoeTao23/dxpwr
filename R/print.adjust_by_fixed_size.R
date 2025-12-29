#' Print method
#'
#' Display a human-readable summary of adjust_by_fixed_size objects containing sample size before and after adjustment.
#'
#' @param x A `adjust_by_fixed_size` object returned by [adjust_for_one_fixed_test()].
#' @return Invisibly returns for `adjust_by_fixed_size` object.
#' @examples
#' n_adj <- adjust_for_one_fixed_test(n_known_test=52, n_unadjusted=77)
#' print(n_adj)
#' n_adj  # same as print(n_adj)
#' @export
print.adjust_fixed_size <- function(x) {

  cat("\033[1mAdjustment Method:\033[0m ", x$method, "\n\n")

  cat("\033[1mCall:\033[0m\n")
  print(x$call)
  cat("\n")

  cat("\033[1mParameters:\033[0m\n")
  cat(sprintf("  %-15s: %d\n", "n_known_test", x$parameter$n_known_test))
  cat(sprintf("  %-15s: %d\n", "n_unadjusted", x$parameter$n_unadjusted))

  cat("\n\033[1mSample Size:\033[0m\n")
  cat(sprintf("  %-6s: n_known_test = %d, n_unknown_test = %d, total = %d\n", "Before", x$sample$n_unadj_known_test, x$sample$n_unadj_unknown_test, x$sample$n_unadjusted))
  cat(sprintf("  %-6s: n_known_test = %d, n_unknown_test = %d, total = %d\n", "After", x$sample$n_adj_known_test, x$sample$n_adj_unknown_test, x$sample$n_adjusted))

  invisible(x)
}
