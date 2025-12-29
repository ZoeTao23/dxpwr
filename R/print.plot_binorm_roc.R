#' Print method
#'
#' Displays a structured summary of `plot_binorm_ROC` objects, including the function call, sample size, and simulation data, followed by the ROC curve plot.
#'
#' @param x A `plot_binorm_roc` object returned by [plot_binorm_ROC()].
#' @return Invisibly returns for `plot_binorm_ROC` object.
#' @examples
#' plot <- plot_binorm_ROC(a=1.19, b=1, n_sample=2000, ci=0.95, seed=1)
#' print(plot)
#' plot  # same as print(plot)
#' @export
print.plot_binorm_roc <- function(x) {
  cat("\033[1mBinormal ROC Plot\033[0m\n")
  cat("\n")

  cat("\033[1mCall:\033[0m\n")
  print(x$call)
  cat("\n")

  cat("\033[1mData Summary:\033[0m\n")
  cat(sprintf("  %-15s: %s\n", "Samples (n)", x$n))
  cat(sprintf("  %-15s: %s\n", "Labels", toString(unique(x$data$labels))))
  cat(sprintf("  %-15s: [%.2f, %.2f]\n", "Scores range", range(x$data$scores)[1], range(x$data$scores)[2]))

  print(x$plot)

  invisible(x)
}
