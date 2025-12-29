#' Sample Size Calculator for Diagnostic Tests
#'
#' Fit a Receiver Operating Characteristic (ROC) curve by assuming normal distribution normally distributed test results for both cases (with condition) and controls (without condition) groups.
#'
#' @param a (numeric): Standardized difference in means of the distributions of test results for patients with and without the condition, (μ_with - μ_without)/σ_with.
#' @param b (numeric): Ratio of the standard deviations of the distributions of test results for patients without versus with the condition, σ_without/σ_with.
#' @param n_sample (Integer): Half sample size for simulation (default: 1000).
#' @param ci_level (numeric):  Confidence level (default: 0.95).
#' @param seed (Integer): Random seed (default: 154).
#' @return A ggplot object showing the ROC curve with the following elements:
#' \itemize{
#'   \item plot (ggplot object) - The ROC curve under binormal model.
#'   \item n (integar) - The total sample size.
#'   \item data (list) - Randomly generated simulation data (scores, labels).
#'   \item call (character) - The original function call.
#' }
#' @importFrom ggplot2 pROC
#' @examples
#' p <- plot_binorm_ROC(a=1.19, b=1, n_sample=2000, ci=0.95, seed=1)
#' print(p)
#' @export
plot_binorm_ROC <- function(a, b, n_sample=1000, ci_level=0.95, seed=154) {

  call <- match.call()

  #  (1) ----- validate inputs
  if (!is.numeric(a) || length(a) != 1) {
    stop("'a' must be a single numeric value.")
  }

  if (!is.numeric(b) || length(b) != 1 || b <= 0) {
    stop("'b' must be a single positive numeric value.")
  }

  if (!is.numeric(ci_level) || length(ci_level) != 1 || ci_level <= 0 || ci_level >= 1) {
    stop("'ci_level' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(n_sample) || length(n_sample) != 1 || n_sample <= 0 || n_sample != round(n_sample)) {
    stop("'n_sample' must be a single positive integer.")
  }

  if (!is.numeric(seed) || length(seed) != 1 ||  seed < 0 || seed != round(seed)) {
    stop("'seed' must be a single positive integer.")
  }

  requireNamespace("ggplot2", quietly = TRUE) ||
    stop("Package 'ggplot2' is required. Please install it first.")
  requireNamespace("pROC", quietly = TRUE) ||
    stop("Package 'pROC' is required. Please install it first.")


  #  (2) ----- calculation modules

  set.seed(seed)

  ## generate binormal samples
  controls <- stats::rnorm(n_sample, mean = 0, sd = 1)
  cases <- stats::rnorm(n_sample, mean = a, sd = b)

  ## construct response labels and scores
  labels <- factor(
    c(rep(0, n_sample), rep(1, n_sample)),
    levels = c(0, 1),
    labels = c("Control", "Case")
  )
  scores <- c(controls, cases)

  ## compute ROC with CI
  roc_obj <- pROC::roc(
    response = labels,
    predictor = scores,
    ci = TRUE,
    quiet = TRUE
  )
  ci_data <- pROC::ci.se(
    roc = roc_obj,
    specificities = seq(0, 1, 0.01)
  )

  ## plot data prepare
  plot_data <- data.frame(
    FPR = 1 - as.numeric(rownames(ci_data)),
    TPR = roc_obj$sensitivities[match(as.numeric(rownames(ci_data)), roc_obj$specificities)],
    Lower = ci_data[, 1],
    Upper = ci_data[, 3]
  )

  ## calculate auc
  auc_val <- round(roc_obj$auc, 3)
  auc_ci <- round(pROC::ci.auc(roc_obj, conf.level = ci_level), 3)


  ## plot
  gg <- ggplot2::ggplot(plot_data) +
    ggplot2::geom_ribbon(
      ggplot2::aes(x=FPR, ymin=Lower, ymax=Upper),
      fill="#293890",
      alpha=0.3
    ) +
    ggplot2::geom_line(
      ggplot2::aes(x=FPR, y=TPR),
      color="#293890",
      linewidth=1
    ) +
    ggplot2::geom_abline(
      slope=1,
      intercept=0,
      linetype="longdash",
      color="#808080",
      linewidth=0.8
    ) +
    ggplot2::annotate(
      "text", x=0.7, y=0.2,
      size=4, fontface="bold",
      label=
        sprintf(
          "AUC = %s\n(%d%% CI: %s-%s)",
          auc_val, ci_level*100,
          auc_ci[1], auc_ci[3]
        )
    ) +
    ggplot2::labs(
      title="ROC Curve under Binormal Model",
      subtitle=sprintf("a = %.2f, b = %.2f, n = %d", a, b, 2*n_sample),
      x="False Positive Rate (1 - Specificity)",
      y="True Positive Rate (Sensitivity)"
    ) +
    ggplot2::theme_bw(base_size=14) +
    ggplot2::theme(
      plot.title =
        ggplot2::element_text(
          hjust=0.5, face="bold", size=14
        ),
      plot.subtitle =
        ggplot2::element_text(
          hjust=0.5, size=12
        ),
      panel.grid.minor =
        ggplot2::element_blank(),
      panel.grid.major =
        ggplot2::element_blank()
    )


  # (3) ----- Structured outputs
  p <-
    structure(
      list(
        plot = gg,
        n = n_sample * 2,
        data = list(
          labels = labels,
          scores = scores
        ),
        call = match.call()
      ),
      class = "plot_binorm_roc"
    )

  return(p)

}
