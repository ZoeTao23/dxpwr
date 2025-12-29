#' Print method
#'
#' Display a human-readable summary of binorm_params objects containing binormal model parameters (a,b) for Receiver Operating Characteristic (ROC) curve fitting.
#'
#' @param x A `binorm_params` object returned by [get_binorm_params()].
#' @param digits Number of significant digits to display (default: 3).
#' @return Invisibly returns for `binorm_params` object.
#' @examples
#' params <- get_binorm_params(AUC = 0.8, SD_Ratio = 1.2)
#' print(params)
#' params  # same as print(params)
#' @export
print.binorm_params <- function(x, digits = 3) {
  cat("\033[1mBinormal ROC Parameters\033[0m\n")
  cat("\n")

  cat("\033[1mCall:\033[0m\n")
  print(x$call)
  cat("\n")

  cat("\033[1mParameters:\033[0m\n")
  cat(sprintf(paste0("  %-", 30, "s = % .", digits, "f\n"),
              "a: (μ_without-μ_with)/σ_with", x[["a"]]))
  cat(sprintf(paste0("  %-", 30, "s = % .", digits, "f\n"),
              "b: σ_without/σ_with", x[["b"]]))

  invisible(x)
}
