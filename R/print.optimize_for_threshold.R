#' Print method
#'
#' Display a human-readable summary of optimize_for_cutoff objects containing the required sample size for deciding a suitable cut-off.
#'
#' @param x A `optimize_for_cutoff` object returned by [optimize_for_threshold()].
#' @return Invisibly returns for `optimize_for_cutoff` object.
#' @examples
#' n <- optimize_for_threshold(spec_min=0.9, sens_at_spec_min=0.6, sens_min=0.5, alpha=0.05, beta=0.2, b=1)
#' print(n)
#' n # same as print(n_adj)
#' @export
print.optimize_for_threshold <- function(x) {

  cat("\033[1mAdjustment Method:\033[0m ", x$method, "\n\n")

  cat("\033[1mCall:\033[0m\n")
  print(x$call)
  cat("\n")

  cat("\033[1mParameters:\033[0m\n")
  cat(sprintf("  %-20s: %.2f\n", "min(spec.)", x$parameter$min_spec))
  cat(sprintf("  %-20s: %.2f\n", "sens. at min(spec.)", x$parameter$sens_at_min_spec))
  cat(sprintf("  %-20s: %.2f\n", "min(sens.)", x$parameter$min_sens))
  cat(sprintf("  %-20s: %.2f\n", "alpha", x$parameter$alpha))
  cat(sprintf("  %-20s: %.2f\n", "power", x$parameter$power))
  cat(sprintf("  %-20s: %.2f\n", "b (binormal params)", x$parameter$b))

  cat("\n\033[1mSample Size:\033[0m\n")
  cat(sprintf("  %-6s: total = %d", "Required", x$sample$n_total))

  invisible(x)
}


