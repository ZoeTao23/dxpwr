#' Sample Size Calculator for Diagnostic Tests
#'
#' Calculate the required sample size for estimating the Area Under the ROC Curve (AUC) of a diagnostic test.
#'
#' @param A (numeric): Expected area under the ROC curve.
#' @param b (numeric): Ratio of the standard deviations of the distributions of test results for patients without versus with the condition, σ_without/σ_with (default: 1).
#' @param alpha (Numeric): Significance level (default: 0.05).
#' @param beta (Numeric): Type II error rate (default: NULL).
#' @param L (Numeric): The desired length of one-half of the (1-α)×100% confidence interval for AUC.
#' @param R (Numeric): Ratio of patients with and without the condition (default: 1).
#' @param dist (Character): The assumption for choosing the variance function (default: "any"):
#'   - `"any"`: Applicable for all distributions.
#'   - `"exp"`: Assume that the unobserved, underlying test results follow an exponential distribution, but the observed test results, either continuous or ordinal, do not necessarily have an exponential distribution
#'   - `"binorm"`: Assume that the unobserved, underlying test results follow a binormal distribution, but the observed test results, either continuous or ordinal, do not necessarily have a binormal distribution
#'   - `"obs_binorm"`: Assume that the observed test results are on a truly continuous scale and follow a binormal distribution or can be transformed to a binormal distribution.
#' @return An object of class "diag_sample_size_html" containing:
#'   - `sample_size` (list): Required sample size.
#'   - `parameters` (list): Input parameters.
#'   - `html_report` (html): Sample size calculation report.
#' @examples
#' lapply(X=c("any", "exp", "binorm", "obs_binorm"), FUN=function(dist) sample_estimate_full_ROC(A=0.8, alpha=0.05, beta=0.2, L=0.1, R=1, dist=dist))
sample_estimate_full_ROC <- function(A, b=1, alpha=0.05, beta=NULL, L, R=1, dist="any") {


  #  (1) ----- validate inputs
  if (!is.numeric(A) || length(A) != 1 || A <= 0 || A >= 1) {
    stop("'A' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(b) || length(b) != 1 || b <= 0) {
    stop("'b' must be a single positive numeric value.")
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value in (0, 1).")
  }

  if (!is.null(beta) && (!is.numeric(beta) || length(beta) != 1 || beta <= 0 || beta >= 1)) {
    stop("'beta' must be a single numeric value in (0, 1), or NULL.")
  }

  if (!is.numeric(L) || length(L) != 1 || L <= 0 || L > 0.5) {
    stop("'L' must be a numeric value in (0, 0.5).")
  }

  if (!is.numeric(R) || length(R) != 1 || R <= 0) {
    stop("'R' must be a single positive numeric value.")
  }

  if (!is.character(dist) || length(dist) != 1 || !(dist %in% c("any", "exp", "binorm", "obs_binorm"))) {
    stop("'dist' must be one of 'any', 'exp', 'binorm' or 'obs_binorm'.")
  }


  #  (2) ----- calculation modules

  ## calculate variance function
  var_function <- switch(
    dist,
    "any" = {
      method <- "Sample size method for estimating AUC by Blume (2009)"
      A * (1 - A)
    },
    "exp" = {
      method <- "Sample size method for estimating AUC by Hanley and McNeil (1982)"
      A / (2 - A) / R + 2 * A^2 / (1 + A) - A^2 * (1 / R + 1)
    },
    "binorm" = {
      method <- "Sample size method for estimating AUC by Obuchowski (1994)"
      a <- get_binorm_params(A, b)$a

      expr1 <- exp(- a^2 / 2 / (1 + b^2))
      expr2 <- 1 + b^2

      f <- expr1 * (2 * pi * expr2)^(-1/2)
      g <- -expr1 * (a*b) * (2 * pi * expr2^3)^(-1/2)

      # (0.0099) * exp(-a^2 / 2) * ((5 * a^2 + 8) + (a^2 + 8) / R) # if b=1
      f^2 * (1 + b^2 / R + a^2 / 2) + g^2 * (b^2 * (1 + R)/ 2 / R)
    },
    "obs_binorm" = {
      method <- "Sample size method for estimating AUC by Obuchowski and McClish (1997)"
      a <- get_binorm_params(A, b)$a

      expr1 <- exp(- a^2 / 2 / (1 + b^2))
      expr2 <- 1 + b^2

      f <- expr1 * (2 * pi * expr2)^(-1/2)
      g <- -expr1 * (a*b) * (2 * pi * expr2^3)^(-1/2)

      # (0.0099) * exp(-a^2 / 2) * ((5 * a^2 + 8) + (a^2 + 8) / R) - 0.0398 * a^2 * exp(-a^2 / 2) # if b = 1
      f^2 * (1 + b^2 / R + a^2 / 2) + g^2 * (b^2 * (1 + R)/ 2 / R) + f*g*a*b

    },
    stop("Invalid 'dist' value. Use 'any', 'exp', 'binorm', or 'obs_binorm'.")
  )

  ## calculate sample size
  N <- get_diag_sample(
    var_function = var_function,
    alpha = alpha,
    beta = beta,
    delta = L,
    test_type = "one_diagnostic"
  )

  ## calculate number of patients with and without the condition
  if (dist == "any" & R < 1) {

    n_with_condition <- ceiling(N*R)
    n_without_condition <- ceiling(N)

  } else {

    n_with_condition <- ceiling(N)
    n_without_condition <- ceiling(N*R)
  }

  n_total <- n_with_condition + n_without_condition


  #  (3) ----- structured outputs
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
        htmltools::tags$strong(style = paste0("color:", ifelse(A > 0.5, "#E74C3C", "#27AE60"), ";"),"area under the ROC curve"),
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
        Parameter = c("Target Accuracy",
                      "Binominal Parameter",
                      "Confidence Level",
                      "Statistical Power",
                      "Margin of Error",
                      "Group Allocation"
        ),
        Symbol    = c("A",
                      "(a, b)",
                      "1 - α",
                      "1 - β",
                      "± L",
                      "R"
        ),
        Value     = c(A,
                      ifelse(dist %in% c("binorm","obs_binorm"), paste0("(", round(a,3), ", ", round(b,3), ")"), "Not specified"),
                      paste0((1-alpha)*100, "%"),
                      ifelse(is.null(beta), "Not specified", paste0((1-beta)*100, "%")),
                      paste0("±", L*100, "%"),
                      R
        ),
        Notes     = c("Desired AUC",
                      "a = (μ_with - μ_without)/σ_with; b = σ_without/σ_with",
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
          "For estimating area under the ROC curve (AUC) with a given precision, the required sample size is calculated using the formula:"
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

        if(dist=="exp") {

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

              paste0("$$ \\hat{V}(\\hat{A})=\\dfrac{A}{R(2-A)}+\\dfrac{2A^2}{1+A}-A^2(\\dfrac{1}{R}+1) $$")
            )
          )

        } else if(dist=="binorm") {

          htmltools::tagList(

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

                paste0("$$ \\hat{V}(\\hat{A})= \\frac{e^{-a^2/2}}{4\\pi(1+b^2)^3} \\times\\{a^2[b^4+(1+b^2)^2]+2(1+b^2)^2+\\frac{a^2b^4+2b^2(1+b^2)^2}{R}\\} $$")
              )
            ),

            htmltools::tags$p(
              "Specifically, when b = 1"
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

                paste0("$$ \\hat{V}(\\hat{A})=0.0099 \\times e^{-a^2/2} \\times(5a^2+8+\\dfrac{a^2+8}{R}) $$")
              )
            )
          )

        } else if(dist=="any") {

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

              paste0("$$ \\hat{V}(\\hat{A})=A(1-A) $$")
            )
          )

        } else if(dist=="obs_binorm") {

          htmltools::tagList(

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

                paste0("$$ \\hat{V}(\\hat{A})=\\frac{e^{-a^2/2}}{4\\pi(1+b^2)^3} \\times[a^2+2(1+b^2)^2+\\frac{a^2b^4+2b^2(1+b^2)^2}{R}] $$")
              )
            ),

            htmltools::tags$p(
              "Specifically, when b = 1"
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

                paste0("$$ \\hat{V}(\\hat{A})=0.0099 \\times e^{-a^2/2} \\times(5a^2+8+\\frac{a^2+8}{R}) $$")
              )
            )

          )

        }


      ),


      ## step 3
      htmltools::tags$div(
        class = "formula-step",

        htmltools::tags$h4(
          style = "color: #2980B9;",
          "Step 3: Calculation Process"
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
                "$$z_{\\alpha/2} = z_{%.3f} = %.3f $$",
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
          "Next, compute the variance function term:"
        ),


        if(dist=="exp") {

          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",
            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
              sprintf(
                paste0("$$ \\hat{V}(\\hat{A}) = \\dfrac{%.3f}{%.3f \\times (2-%.3f)}+\\dfrac{2\\times%.3f^2}{1+%.3f}-%.3f^2\\times(\\dfrac{1}{%.3f}+1)= %.3f $$"),
                A, R, A, A, A, A, R, var_function
              )
            )
          )


        } else if(dist=="binorm") {

          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",
            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
              sprintf(
                paste0("$$ \\hat{V}(\\hat{A}) = %.4f \\times e^{-%.3f^2/2} \\times (%.3f \\times%.3f^2+ %.3f +\\dfrac{%.3f \\times %.3f^2+%.3f}{%.3f})= %.3f $$"),
                1/32/pi, a, b^4+(1+b^2)^2, a, 2*(1+b^2)^2, b^4, a, 2*b^2*(1+b^2)^2, R, var_function
              )
            )
          )

        } else if(dist=="any") {

          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",
            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
              sprintf(
                paste0("$$ \\hat{V}(\\hat{A}) =%.3f\\times(1-%.3f) = %.3f $$"),
                A, A, var_function
              )
            )
          )

        } else if(dist=="obs_binorm") {

          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",
            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
              sprintf(
                paste0("$$ \\hat{V}(\\hat{A}) = %.4f \\times e^{-%.3f^2/2} \\times (%.3f^2+ %.3f +\\dfrac{%.3f \\times %.3f^2+%.3f}{%.3f})= %.3f  $$"),
                1/32/pi, a, a, 2*(1+b^2)^2, b^4, a, 2*b^2*(1+b^2)^2, R, var_function
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
                  "$$n = \\dfrac{[(%.3f)\\times \\sqrt{%.3f}]^2}{(%.3f)^2} ≈ %d$$",
                  qnorm(alpha/2),
                  var_function,
                  L,
                  N
                )
              )

            )

          )
        } else {
          htmltools::tags$div(
            style = "text-align: center; margin: 15px 0;",

            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",

              htmltools::HTML(
                sprintf(
                  "$$n = \\dfrac{\\{[(%.3f)+(%.3f)]\\times \\sqrt{%.3f}\\}^2}{(%.3f)^2} ≈ %d$$",
                  qnorm(alpha/2),
                  qnorm(beta),
                  var_function,
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
        "Sample size calculations for estimating area under the ROC curve"
      ),

      ## 2.Assumptions
      if(dist=="exp") {

        htmltools::tags$p(
          style = "margin: 15px 0px; text-indent: -1em; padding-left: 1.5em;",
          htmltools::tags$strong("2. Assumptions:"),
          htmltools::tags$br(),
          "Assume that the unobserved, underlying test results follow an exponential distribution, but the observed test results, either continuous or ordinal, do not necessarily have an exponential distribution"

        )

      } else if(dist=="binorm") {

        htmltools::tags$p(
          style = "margin: 15px 0px; text-indent: -1em; padding-left: 1.5em;",
          htmltools::tags$strong("2. Assumptions:"),
          htmltools::tags$br(),
          "Assume that the unobserved, underlying test results follow a binormal distribution, but the observed test results, either continuous or ordinal, do not necessarily have a binormal distribution"

        )

      } else if(dist=="any") {

        htmltools::tags$p(
          style = "margin: 15px 0px; text-indent: -1em; padding-left: 1.5em;",
          htmltools::tags$strong("2. Assumptions:"),
          htmltools::tags$br(),
          "The sample size that would be applicable for the distribution that would yield the smallest sample size and the smallest sample size that would be applicable for all distributions"
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
        if(dist=="exp") {

          htmltools::tagList(
            "[3] Hanley JA, McNeil BJ. The Meaning and Use of the Area under a Receiver Operating Characteristic (ROC) curve. ",
            htmltools::tags$em("Radiology."),
            " 1982;143(1):29-36. ",
            htmltools::tags$a(
              href = "https://doi.org/10.1148/radiology.143.1.7063747",
              target = "_blank",
              style = "color: #1a5fb4; text-decoration: none; border-bottom: 1px solid #c1dbf7;
               transition: border-color 0.3s, color 0.3s;",
              onmouseover = "this.style.borderBottomColor='#1a5fb4'; this.style.color='#0d4e9b';",
              onmouseout = "this.style.borderBottomColor='#c1dbf7'; this.style.color='#1a5fb4';",
              "https://doi.org/10.1148/radiology.143.1.7063747"
            )
          )

        } else if(dist=="binorm") {

          htmltools::tagList(
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

        } else if(dist=="any") {

          htmltools::tagList(
            "[3] Blume JD. Bounding Sample Size Projections for the Area under a ROC Curve. ",
            htmltools::tags$em("Journal of Statistical Planning and Inference."),
            " 2009;139(3):711-721. ",
            htmltools::tags$a(
              href = "doi:https://doi.org/10.1016/j.jspi.2007.09.015",
              target = "_blank",
              style = "color: #1a5fb4; text-decoration: none; border-bottom: 1px solid #c1dbf7;
               transition: border-color 0.3s, color 0.3s;",
              onmouseover = "this.style.borderBottomColor='#1a5fb4'; this.style.color='#0d4e9b';",
              onmouseout = "this.style.borderBottomColor='#c1dbf7'; this.style.color='#1a5fb4';",
              "doi:https://doi.org/10.1016/j.jspi.2007.09.015"
            )
          )

        } else if(dist=="obs_binorm") {

          htmltools::tagList(
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
            ),

            htmltools::tags$br(),
            "[4] Obuchowski NA, McClish DK. Sample Size Determination for Diagnostic Accuracy Studies Involving Binormal ROC Curve Indices. ",
            htmltools::tags$em("Statistics in Medicine."),
            " 1997;16(13):1529-1542. ",
            htmltools::tags$a(
              href = "doi:https://doi.org/10.1002/(sici)1097-0258(19970715)16:13%3C1529::aid-sim565%3E3.0.co;2-h",
              target = "_blank",
              style = "color: #1a5fb4; text-decoration: none; border-bottom: 1px solid #c1dbf7;
               transition: border-color 0.3s, color 0.3s;",
              onmouseover = "this.style.borderBottomColor='#1a5fb4'; this.style.color='#0d4e9b';",
              onmouseout = "this.style.borderBottomColor='#c1dbf7'; this.style.color='#1a5fb4';",
              "doi:https://doi.org/10.1002/(sici)1097-0258(19970715)16:13%3C1529::aid-sim565%3E3.0.co;2-h"
            )
          )

        }

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

