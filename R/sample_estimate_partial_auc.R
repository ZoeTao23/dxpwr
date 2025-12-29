#' Sample Size Calculator for Diagnostic Tests
#'
#'Calculate the required sample size for estimating the partial Area Under the ROC Curve (pAUC) of a diagnostic test.
#'
#' @param A (numeric): Expected area under the ROC curve.
#' @param e1 (numeric): Lower bound of the false positive rate (FPR) range.
#' @param e2 (numeric): Upper bound of the false positive rate (FPR) range.
#' @param b (numeric): Ratio of the standard deviations of the distributions of test results for patients without versus with the condition, σ_without/σ_with (default: 1).
#' @param alpha (numeric): Significance level (default: 0.05).
#' @param beta (numeric): Type II error rate (default: NULL).
#' @param L (numeric): The desired length of one-half of the (1-α)×100% confidence interval for partial AUC.
#' @param R (Numeric): Ratio of patients with and without the condition (default: 1).
#' @param dist (character): The assumption for choosing the variance function (default: "binorm"):
#'   - `"binorm"`: Assume that the unobserved, underlying test results follow a binormal distribution, but the observed test results, either continuous or ordinal, do not necessarily have a binormal distribution.
#'   - `"obs_binorm"`: The observed test results follow a binormal distribution.
#' @return An object of class "diag_sample_size_html" containing:
#'   - `sample_size` (list): Required sample size.
#'   - `parameters` (list): Input parameters.
#'   - `html_report` (html): Sample size calculation report.
#' @import stats
#' @examples
#' n <- sample_estimate_partial_auc(A=0.8, e1=0, e2=0.1, b=1, alpha=0.05, beta=NULL, L=0.1, R=1, dist="binorm")
#' print(n)
#' @export
sample_estimate_partial_auc <- function(A, b=1, e1, e2, alpha=0.05, beta=NULL, L, R=1, dist="binorm") {


  #  (1) ----- validate inputs

  if (!is.numeric(A) || length(A) != 1 || A <= 0 || A >= 1) {
    stop("'A' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(b) || length(b) != 1 || b <= 0) {
    stop("'b' must be a single positive numeric value.")
  }

  if (!is.numeric(e1) || length(e1) != 1 || e1 < 0 || e1 >= 1) {
    stop("'e1' must be a single numeric value in [0, 1).")
  }

  if (!is.numeric(e2) || length(e2) != 1 || e2 <= 0 || e2 > 1) {
    stop("'e2' must be a single numeric value in (0, 1].")
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric in (0, 1).")
  }

  if (!is.null(beta) && (!is.numeric(beta) || length(beta) != 1 || beta <= 0 || beta >= 1)) {
    stop("'beta' must be a single numeric value in (0, 1), or NULL.")
  }

  if (!is.numeric(L) || length(L) != 1 || L <= 0 || L > 0.5) {
    stop("'L' must be a single numeric value in (0, 0.5).")
  }

  if (!is.numeric(R) || length(R) != 1 || R <= 0) {
    stop("'R' must be a single positive numeric value.")
  }

  if (!is.character(dist) || length(dist) != 1 || !(dist %in% c("binorm", "obs_binorm"))) {
    stop("'dist' must be either 'binorm' or 'obs_binorm'.")
  }


  #  (2) ----- calculation modules


  ## calculate parameters a, e', and e''
  a <- get_binorm_params(A, b)$a
  e_prime <- (qnorm(c(e1,e2)) + a * b * (1+b^2)^(-1)) * sqrt(1 + b^2)
  e_double_prime <- (e_prime)^2 / 2


  ## calculate variance function
  expr1 <- exp(-a^2 / 2 / (1 + b^2))
  expr2 <- 1 + b^2
  expr3 <- pnorm(e_prime[2]) - pnorm(e_prime[1])
  expr4 <- exp(-e_double_prime[1]) - exp(-e_double_prime[2])

  f <- expr1 * (2 * pi * expr2)^(-1 / 2) * expr3
  g <- expr1 * (2 * pi * expr2)^(-1) * expr4 - a * b * expr1 * (2 * pi * expr2^3)^(-1 / 2) * expr3

  var_function <- switch(
    dist,
    "binorm" = f^2 * (1 + b^2 / R + a^2 / 2) + g^2 * (b^2 * (1 + R) / (2 * R)),
    "obs_binorm" = f^2 * (1 + b^2 / R + a^2 / 2) + g^2 * (b^2 * (1 + R) / (2 * R)) + f * g * a * b,
    stop("Invalid 'dist' value. Use 'binorm' or 'obs_binorm'.")
  )

  ## calculate delta
  delta <- abs((e1-e2)*L)

  ## calculate sample size
  N <- get_diag_sample(
    var_function = var_function,
    alpha = alpha,
    beta = beta,
    delta = delta,
    test_type = "one_diagnostic"
  )

  ## calculate number of patients with and without condition
  n_with_condition <- ceiling(N)
  n_without_condition <- ceiling(N*R)
  n_total <- n_with_condition + n_without_condition


  #  (3) ----- structured outputs

  ## method reference
  html_report <- htmltools::tags$div(
    class = "clinical-report",
    style = "font-family: Arial; max-width: 800px; margin: auto;",


    # @headline
    htmltools::tags$h2(
      style = "color: #2E86C1; border-bottom: 2px solid #3498DB; padding-bottom: 10px;",
      "Diagnostic Accuracy Study Sample Size Calculation"
    ),


    # @study summary
    htmltools::tags$div(
      style = "background-color:#F8F9F9; padding: 15px; border-radius: 5px; margin-bottom: 20px; line-height: 1;",
      htmltools::tags$h3(style="margin-top: 0;", "🖥️ Study Summary"),

      ## 1.Objective
      htmltools::tags$p(
        htmltools::tags$strong("Objective:"),
        "To evaluate the ",
        htmltools::tags$strong(style = paste0("color:", ifelse(A > 0.5, "#E74C3C", "#27AE60"), ";"),"partial area under the ROC curve"),
        "of a new diagnostic test"
      ),

      ## 2.Type
      htmltools::tags$p(
        htmltools::tags$strong("Type:"),
        "Clinical performance evaluation"
      ),

      ## 3.Phase
      htmltools::tags$p(
        htmltools::tags$strong("Phase:"),
        "Phase I – exploration, Phase II – development and Phase III – validation"
      )

    ),


    # @key parameters
    htmltools::tags$h3("📊 Key Parameters"),
    DT::datatable(
      data.frame(
        Parameter = c("Area under the ROC Curve",
                      "Binominal Parameter",
                      "FPR Range",
                      "Confidence Level",
                      "Statistical Power",
                      "Margin of Error",
                      "Group Allocation"
        ),
        Symbol    = c("A",
                      "(a, b)",
                      "(e1, e2)",
                      "1 - α",
                      "1 - β",
                      "± L",
                      "R"
        ),
        Value     = c(A,
                      paste0("(", round(a,3), ", ", round(b,3), ")"),
                      paste0("(", round(e1,3), ", ", round(e2,3), ")"),
                      paste0((1-alpha)*100, "%"),
                      ifelse(is.null(beta), "Not specified", paste0((1-beta)*100, "%")),
                      paste0("±", L*100, "%"),
                      R
        ),
        Notes     = c("Desired AUC, to determine the shape of ROC",
                      "a = (μ_with - μ_without)/σ_with; b = σ_without/σ_with, to determine the shape of ROC",
                      "e1 = Lower bound of FPR, e2 = Upper bound of FPR",
                      "1 - Type I error rate",
                      "1 - Type II error rate",
                      "Half-width of CI",
                      "Ratio of patients with and without the condition"
        )
      ),

      options = list(dom = 't', pageLength = 8),
      rownames = FALSE,
      class = "stripe hover"
    ),


    # @formula derivation
    htmltools::tags$h3("🔬 Formula Derivation"),
    htmltools::tags$div(
      class = "formula-container",
      style = "background-color: #f8f9fa; border-radius: 8px; padding: 20px;",

      htmltools::tags$div(
        class = "formula-step",

        ## step 1
        htmltools::tags$h4(
          style = "color: #2980B9; margin-top: 0;",
          "Step 1: Sample Size Formula"
        ),

        htmltools::tags$p(
          "For estimating partial area under the ROC curve (AUC) with a given precision, the required sample size is calculated using the formula:"
        ),

        if (is.null(beta)) {
          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",

            htmltools::tags$p(
              style = paste(
                "font-family: 'Cambria Math', 'Latin Modern Math', serif;",
                # "font-size: 1.2em;",
                "padding: 5px;",
                "background-color: white;",
                "border-radius: 5px;",
                "box-shadow: 0 2px 5px rgba(0,0,0,0.05);",
                "display: inline-block;"
              ),

              "$$ n = \\dfrac{[Z_{\\alpha/2} \\sqrt{V(\\hat{A})}]^2}{L^2} $$"
            )
          )

        } else {
          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",

            htmltools::tags$p(
              style = paste(
                "font-family: 'Cambria Math', 'Latin Modern Math', serif;",
                # "font-size: 1.2em;",
                "padding: 5px;",
                "background-color: white;",
                "border-radius: 5px;",
                "box-shadow: 0 2px 5px rgba(0,0,0,0.05);",
                "display: inline-block;"
              ),

              "$$ n = \\dfrac{[(Z_{\\alpha/2}+Z_{\\beta}) \\sqrt{V(\\hat{A})}]^2}{L^2} $$"
            )
          )

        }

      ),

      ## step 2
      htmltools::tags$div(
        class = "formula-step",

        htmltools::tags$h4(
          style = "color: #2980B9; margin-top: 0;",
          "Step 2: Variance Function"
        ),

        htmltools::tags$p(
          "The variance function is calculated through:"
        ),

        if(dist=="binorm") {

          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",

            htmltools::tags$p(
              style = paste(
                "font-family: 'Cambria Math', 'Latin Modern Math', serif;",
                # "font-size: 1.2em;",
                "padding: 5px;",
                "background-color: white;",
                "border-radius: 5px;",
                "box-shadow: 0 2px 5px rgba(0,0,0,0.05);",
                "display: inline-block;"
              ),

              paste0("$$ \\hat{V}(\\hat{A}_{e_{1} \\leq FPR \\leq e_{2}})= f^2 (1+\\frac{b^2}{R}+\\frac{a^2}{2})+g^2(b^2\\frac{1+R}{2R}) $$")
            )
          )

          } else if(dist=="obs_binorm") {

          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",

            htmltools::tags$p(
              style = paste(
                "font-family: 'Cambria Math', 'Latin Modern Math', serif;",
                # "font-size: 1.2em;",
                "padding: 5px;",
                "background-color: white;",
                "border-radius: 5px;",
                "box-shadow: 0 2px 5px rgba(0,0,0,0.05);",
                "display: inline-block;"
              ),

              paste0("$$ \\hat{V}(\\hat{A}_{e_{1} \\leq FPR \\leq e_{2}})= f^2 (1+\\frac{b^2}{R}+\\frac{a^2}{2})+g^2(b^2\\frac{1+R}{2R})+fgab $$")
            )
          )

        },

        htmltools::tags$p(
          "Where f is a function of (a, b):"
        ),

        htmltools::tags$div(
          style = "text-align: center; margin: 10px 0;",

          htmltools::tags$p(
            style = paste(
              "font-family: 'Cambria Math', 'Latin Modern Math', serif;",
              # "font-size: 1.2em;",
              "padding: 5px;",
              "background-color: white;",
              "border-radius: 5px;",
              "box-shadow: 0 2px 5px rgba(0,0,0,0.05);",
              "display: inline-block;"
            ),

            paste0("$$ f = e^{-\\frac{a^2}{2(1+b^2)}} \\times \\frac{\\Phi(e'_{2}) -\\Phi(e'_{1})}{\\sqrt{2\\pi (1+b^2)}} $$")
          )
        ),

        htmltools::tags$p(
          "And similarly for g:"
        ),

        htmltools::tags$div(
          style = "text-align: center; margin: 10px 0;",

          htmltools::tags$p(
            style = paste(
              "font-family: 'Cambria Math', 'Latin Modern Math', serif;",
              # "font-size: 1.2em;",
              "padding: 5px;",
              "background-color: white;",
              "border-radius: 5px;",
              "box-shadow: 0 2px 5px rgba(0,0,0,0.05);",
              "display: inline-block;"
            ),

            paste0("$$ g = e^{-\\frac{a^2}{2(1+b^2)}} \\times [\\frac{e^{-(e'_{2})^2/2} -e^{-(e'_{1})^2/2}}{\\sqrt{2\\pi (1+b^2)}} - \\frac{ab[\\Phi(e'_{2})-\\Phi(e'_{1})]}{\\sqrt{2\\pi(1+b^2)^3}}] $$")
          )
        ),


        htmltools::tags$p(
          "Where"
        ),

        htmltools::tags$div(
          style = "text-align: center; margin: 10px 0;",
          htmltools::tags$p(
            style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
            "$$ e'_{i} = \\sqrt{1+b^2}[\\Phi^{-1}(e_{i})+\\frac{ab}{1+b^2}] $$",
          )
        )

      ),


      ## step 3
      htmltools::tags$div(
        class = "formula-step",

        htmltools::tags$h4(
          style = "color: #2980B9;",
          "Step 3: Transformed Margin of Error"
        ),

        htmltools::tags$p(
          "In terms of the partial area under the ROC curve"
        ),

        htmltools::tags$p(
          style = "text-align: center; margin: 10px 0; font-family: 'Cambria Math', serif; font-size: 1.1em;",
          paste0("$$ L(\\hat{A}_{e_{1} \\leq FPR \\leq e_{2}}) = L \\times (e_{2}-e_{1}) $$")
        )

      ),

      ## step 4
      htmltools::tags$div(
        class = "formula-step",

        htmltools::tags$h4(
          style = "color: #2980B9;",
          "Step 4: Calculation Process"
        ),

        htmltools::tags$p(
          "First, calculate the Z-score for the given confidence level:"
        ),

        if(is.null(beta)) {
          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",
            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
              sprintf(
                "$$z_{\\alpha/2} = z_{%.2f} = %.3f $$",
                alpha/2, qnorm(alpha/2)
              )
            )
          )
        } else {
          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",
            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
              sprintf(
                "$$z_{\\alpha/2} = z_{%.3f} = %.3f $$",
                alpha/2, qnorm(alpha/2)
              ),
              sprintf(
                "$$z_{\\beta} = z_{%.2f} = %.3f $$",
                beta, qnorm(beta)
              )
            )
          )

        },


        htmltools::tags$p(
          "Next, compute the value of f and g:"
        ),

        htmltools::tags$div(
          style = "text-align: center; margin: 10px 0;",
          htmltools::tags$p(
            style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
            sprintf(
              paste0("$$ e'_{1} = \\sqrt{1+%.3f^2}[\\Phi^{-1}(%.3f)+\\frac{%.3f \\times %.3f}{1+%.3f^2}] = %.3f $$"),
              b, e1, a, b, b, e_prime[1]
            )
          )
        ),

        htmltools::tags$div(
          style = "text-align: center; margin: 10px 0;",
          htmltools::tags$p(
            style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
            sprintf(
              paste0("$$ e'_{2} = \\sqrt{1+%.3f^2}[\\Phi^{-1}(%.3f)+\\frac{%.3f \\times %.3f}{1+%.3f^2}] = %.3f $$"),
              b, e2, a, b, b, e_prime[2]
            )
          )
        ),

        htmltools::tags$p(
          "Plugging in gives f:"
        ),

        htmltools::tags$div(
          style = "text-align: center; margin: 10px 0;",
          htmltools::tags$p(
            style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
            sprintf(
              paste0("$$ f = e^{-\\frac{%.3f^2}{2(1+%.3f^2)}} \\times \\frac{\\Phi(%.3f) -\\Phi(%.3f)}{\\sqrt{2\\pi (1+%.3f^2)}} = %.3f $$"),
              a, b, e_prime[2], e_prime[1], b, f
            )
          )
        ),

        htmltools::tags$p(
          "Similarly, obtain g:"
        ),

        htmltools::tags$div(
          style = "text-align: center; margin: 10px 0;",
          htmltools::tags$p(
            style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
            sprintf(
              paste0("$$ g = e^{-\\frac{%.3f^2}{2\\times(1+%.3f^2)}} \\times [\\frac{e^{-(%.3f)^2/2} -e^{-(%.3f)^2/2}}{\\sqrt{2\\pi (1+%.3f^2)}} - \\frac{%.3f \\times %.3f \\times [\\Phi(%.3f)-\\Phi(%.3f)]}{\\sqrt{2\\pi(1+%.3f^2)^3}} = %.3f $$"),
              a, b, e_prime[2], e_prime[1], b, a, b, e_prime[2], e_prime[1], b, g
            )
          )
        ),


        htmltools::tags$p(
          "Next, compute the variance function term:"
        ),

        if(dist=="binorm") {

          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",
            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
              sprintf(
                paste0("$$ \\hat{V}(\\hat{A}_{e_{1} \\leq FPR \\leq e_{2}})= %.3f^2 \\times (1+\\frac{%.3f^2}{%.3f}+\\frac{%.3f^2}{2})+(%.3f)^2 \\times (%.3f^2 \\times \\frac{1+%.3f}{2 \\times %.3f}) = %.3f $$"),
                f, b, R, a, g, b, R, R, var_function
              )
            )
          )

        } else if(dist=="obs_binorm") {

          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",
            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
              sprintf(
                paste0("$$\\hat{V}(\\hat{A}_{e_{1} \\leq FPR \\leq e_{2}})= %.3f^2 \\times (1+\\frac{%.3f^2}{%.3f}+\\frac{%.3f^2}{2})+(%.3f)^2 \\times (%.3f^2\\frac{1+%.3f}{2 \\times %.3f})+%.3f \\times %.3f \\times %.3f \\times %.3f = %.3f $$"),
                f, b, R, a, g, b, R, R, f, g, a, b, var_function
              )
            )
          )

        },


        htmltools::tags$p(
          "Now apply the full formula:"
        ),

        if(is.null(beta)) {
          htmltools::tags$div(
            style = "text-align: center; margin: 15px 0;",

            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",

              htmltools::HTML(
                sprintf(
                  "$$n = \\dfrac{[(%.3f)\\times \\sqrt{%.3f}]^2}{[(%.3f-%.3f) \\times %.3f]^2} ≈ %d$$",
                  qnorm(alpha/2),
                  var_function,
                  e2,
                  e1,
                  L,
                  N
                )
              )

            )

          )
        } else{
          htmltools::tags$div(
            style = "text-align: center; margin: 15px 0;",

            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",

              htmltools::HTML(
                sprintf(
                  "$$n = \\dfrac{\\{[(%.3f)+(%.3f)]\\times \\sqrt{%.3f}\\}^2}{[(%.3f+%.3f) \\times %.3f]^2} ≈ %d$$",
                  qnorm(alpha/2),
                  qnorm(beta),
                  var_function,
                  e2,
                  e1,
                  L,
                  N
                )
              )

            )

          )
        },

        htmltools::tags$p(
          "Thus,"
        ),

        htmltools::tags$div(
          style = "text-align: center; margin: 10px 0;",
          htmltools::tags$p(
            style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
            sprintf(
              paste0("$$ N_{+} = n = %d, N_{-} = R \\times n = %d$$"),
              n_with_condition, n_without_condition
            )
          )
        )


      ),


      htmltools::tags$script(
        src = "https://polyfill.io/v3/polyfill.min.js?features=es6"
      ),
      htmltools::tags$script(
        id = "MathJax-script",
        async = TRUE,
        src = "https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js"
      ),
      htmltools::tags$script(
        htmltools::HTML(
          "MathJax = {
            tex: {
              inlineMath: [['$', '$'], ['\\\\(', '\\\\)']]
            },
            svg: {
              fontCache: 'global'
            }
          };"
        )
      )

    ),


    # @sample size requirement
    htmltools::tags$h3("🧮 Sample Size Requirement"),
    htmltools::tags$div(
      style = "background-color: #EBF5FB; padding: 15px; border-left: 4px solid #3498DB;",

      htmltools::tags$h4(style = "margin-top: 0;", "Minimum Participants"),

      htmltools::tags$p(
        style = "font-size: 24px; font-weight: bold; color: #2E86C1;",
        ceiling(n_total)
      ),

      htmltools::tags$h4(style = "margin-top: 0;", "Participants with / without condition"),

      htmltools::tags$p(
        style = "font-size: 24px; font-weight: bold; color: #2E86C1;",
        paste0(n_with_condition," / ", n_without_condition)
      ),

      htmltools::tags$p(
        paste0("This provides ", (1-alpha)*100, "% confidence that the estimated AUC"),
        paste0(" will be within ±", L*100, "% of the true value"),
        ifelse(!is.null(beta), paste0(" with ", (1-beta)*100, "% power"), "")
      )
    ),


    # @Methodology
    htmltools::tags$h3("📚 Methodological References"),
    htmltools::tags$div(
      style = "padding: 0 20px;",

      ## 1.Method
      htmltools::tags$p(
        style = "margin: 15px 0px; text-indent: -1em; padding-left: 1.5em;",
        htmltools::tags$strong("1. Method:"),
        htmltools::tags$br(),
        "Sample size calculations for estimating partial area under the ROC curve"
      ),

      ## 2.Assumptions
      if(dist=="binorm") {

        htmltools::tags$p(
          style = "margin: 15px 0px; text-indent: -1em; padding-left: 1.5em;",
          htmltools::tags$strong("2. Assumptions:"),
          htmltools::tags$br(),
          "Assume that the unobserved, underlying test results follow a binormal distribution, but the observed test results, either continuous or ordinal, do not necessarily have a binormal distribution"

        )

      } else if(dist=="obs_binorm") {


        htmltools::tags$p(
          style = "margin: 15px 0px; text-indent: -1em; padding-left: 1.5em;",
          htmltools::tags$strong("2. Assumptions:"),
          htmltools::tags$br(),
          "Assume that the observed test results are on a truly continuous scale and follow a binormal distribution or can be transformed to a binormal distribution",

        )

      },


      ## 3.Limitations
      htmltools::tags$p(
        style = "margin: 15px 0px; text-indent: -1em; padding-left: 1.5em;",
        htmltools::tags$strong("3. Limitations:"),
        htmltools::tags$br(),
        "Valid on the specific assumptions"
      ),

      ## 4.Reference:
      htmltools::tags$p(
        style = "margin: 15px 0px; line-height: 1.2;text-indent: -1em; padding-left: 1.5em;",
        htmltools::tags$strong("4. Reference:"),

        htmltools::tags$br(),
        "[1] Flahault A, Cadilhac M, Thomas G. Sample Size Calculation Should Be Performed for Design Accuracy in Diagnostic Test Studies. ",
        htmltools::tags$em("Journal of Clinical Epidemiology."),
        " 2005;58(8):859-862. ",
        htmltools::tags$a(
          href = "doi:https://doi.org/10.1016/j.jclinepi.2004.12.009",
          target = "_blank",
          style = "color: #1a5fb4; text-decoration: none; border-bottom: 1px solid #c1dbf7;
             transition: border-color 0.3s, color 0.3s;",
          onmouseover = "this.style.borderBottomColor='#1a5fb4'; this.style.color='#0d4e9b';",
          onmouseout = "this.style.borderBottomColor='#c1dbf7'; this.style.color='#1a5fb4';",
          "doi:https://doi.org/10.1016/j.jclinepi.2004.12.009"
        ),


        htmltools::tags$br(),
        "[2] Neuhaus J, McCulloch C. Generalized Linear Models. ",
        htmltools::tags$em("Wiley Interdisciplinary Reviews: Computational Statistics."),
        " 2011;3(5):407-413. ",
        htmltools::tags$a(
          href = "doi:https://doi.org/10.1002/wics.175",
          target = "_blank",
          style = "color: #1a5fb4; text-decoration: none; border-bottom: 1px solid #c1dbf7;
             transition: border-color 0.3s, color 0.3s;",
          onmouseover = "this.style.borderBottomColor='#1a5fb4'; this.style.color='#0d4e9b';",
          onmouseout = "this.style.borderBottomColor='#c1dbf7'; this.style.color='#1a5fb4';",
          "doi:https://doi.org/10.1002/wics.175"
        ),

        htmltools::tags$br(),
        "[3] Obuchowski NA. Computing Sample Size for Receiver Operating Characteristic Studies. ",
        htmltools::tags$em("Investigative Radiology."),
        " 1994;29(2):238-243. ",
        htmltools::tags$a(
          href = "doi:https://doi.org/10.1097/00004424-199402000-00020",
          target = "_blank",
          style = "color: #1a5fb4; text-decoration: none; border-bottom: 1px solid #c1dbf7;
           transition: border-color 0.3s, color 0.3s;",
          onmouseover = "this.style.borderBottomColor='#1a5fb4'; this.style.color='#0d4e9b';",
          onmouseout = "this.style.borderBottomColor='#c1dbf7'; this.style.color='#1a5fb4';",
          "doi:https://doi.org/10.1097/00004424-199402000-00020"
        )


      )

    ),


    # @Compliance statement
    htmltools::tags$h3("🛡️ Considerations"),
    htmltools::tags$div(
      style = "background-color: #F8F9F9; padding: 1px 20px; border-radius: 5px;",

      htmltools::tags$p(
        style = "line-height: 1;",

        ## 1. Input
        htmltools::tags$p(
          style = "text-indent: -1em; padding-left: 1.5em;",
          htmltools::tags$strong("1. Input Parameters Must Be Evidence-Based"),
          htmltools::tags$br(),
          "All input parameters (e.g., AUC, alpha, beta) should be grounded in real-world data, prior studies, or authoritative literature. We recommend evaluating the impact of ",
          htmltools::tags$u("±20% variation"),
          " in key parameters to ensure robustness of the study design."
        ),

        ## 2. Results
        htmltools::tags$p(
          style = "text-indent: -1em; padding-left: 1.5em;",
          htmltools::tags$strong("2. Results Reflect Minimum Sample Size Requirements"),
          htmltools::tags$br(),
          "The calculated sample size represents ",
          htmltools::tags$u("the minimum number of participants required "),
          "to evaluate the primary diagnostic performance metrics (e.g., AUC and partial AUC) under the specified assumptions."
        ),

        ## 3. Adjustment
        htmltools::tags$p(
          style = "text-indent: -1em; padding-left: 1.5em;",
          htmltools::tags$strong("3. Allow Additional Sample Size Buffer"),
          htmltools::tags$br(),
          "To account for potential issues during the study, it is recommended to increase the calculated sample size ",
          htmltools::tags$u("by 10-15% "),
          "to cover: ",
          htmltools::tags$br(),
          "\u00A0 \u00A0 \u00A0 \u2022 \u00A0 Participant dropout or non-compliance",
          htmltools::tags$br(),
          "\u00A0 \u00A0 \u00A0 \u2022 \u00A0 Protocol deviations",
          htmltools::tags$br(),
          "\u00A0 \u00A0 \u00A0 \u2022 \u00A0 Missing or non-analyzable data"
        ),

        ## 4. Compliance
        htmltools::tags$p(
          style = "text-indent: -1em; padding-left: 1.5em;",
          htmltools::tags$strong("4. Regulatory Compliance"),
          htmltools::tags$br(),
          "For studies intended for regulatory submissions (e.g., FDA, CE Mark, NMPA), ",
          htmltools::tags$u("all assumptions used"),
          " in the sample size calculation must be clearly justified in the Statistical Analysis Plan (SAP)."
        ),

        ## 5. Guidelines
        htmltools::tags$p(
          style = "text-indent: -1em; padding-left: 1.5em;",
          htmltools::tags$strong("5. Align with International Guidelines and Ethical Constraints"),
          htmltools::tags$br(),
          "When planning the study, consider the following:",
          htmltools::tags$br(),
          "\u00A0 \u00A0 \u00A0 \u2022 \u00A0 Minimum statistical power requirements",
          htmltools::tags$br(),
          "\u00A0 \u00A0 \u00A0 \u2022 \u00A0 Maximum sample size permitted by the ethics committee",
          htmltools::tags$br(),
          "\u00A0 \u00A0 \u00A0 \u2022 \u00A0 Relevant international standards, such as:",
          htmltools::tags$br(),
          "\u00A0 \u00A0 \u00A0 \u00A0 \u00A0 \u2022 ", htmltools::tags$u("CLSI EP24-A2"),
          ": Recommends ≥100 reference standard–confirmed positive cases and ≥100 validated disease-free controls",
          htmltools::tags$br(),
          "\u00A0 \u00A0 \u00A0 \u00A0 \u00A0 \u2022 ", htmltools::tags$u("ICH E9, Section 3.5"),
          ": Suggests ≥200 confirmed endpoint events for Phase III clinical trials"
        ),

      ),

      htmltools::tags$p(
        style = "font-style: italic; color: #7D6608; background-color: #FCF3CF; padding: 8px; border-left: 4px solid #F1C40F;",
        "❗  Important: This sample size estimate is for preliminary guidance only. Final determination must be confirmed by certified biostatisticians and approved by the ethics committee and regulatory authorities."        )

    ),


    # @download
    htmltools::tags$div(
      style = "text-align: center; margin-top: 30px;",

      htmltools::tags$button(
        id = "download_report",
        class = "btn btn-primary",
        style = "background-color: #2E86C1; color: white; border: none; padding: 10px 20px; border-radius: 5px;",
        "📥 Download Full Report (PDF)"
      )
    ),


    # 7. JavaScript interaction
    htmltools::tags$script(htmltools::HTML(
      "$(document).ready(function() {
        $('#download_report').click(function() {
          alert('PDF export functionality would be implemented here in production');
        });
      });"
    ))
  )

  # show in Viewer
  if (interactive()) {
    htmltools::html_print(html_report)
  }

  # return the structured object
  invisible(
    structure(
      list(
        sample_size = list(n_total = n_total, n_with_condition = n_with_condition, n_without_condition = n_without_condition),
        parameters = list(A = A, a = ifelse(dist %in% c("binorm","obs_binorm"), a, NA), b = ifelse(dist %in% c("binorm","obs_binorm"), b, NA), alpha = alpha, beta = beta, L = L, R = R),
        html_report = html_report
      ),
      class = "diag_sample_size_html"
    )
  )
}
