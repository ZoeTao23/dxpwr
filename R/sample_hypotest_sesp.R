#' Sample Size Calculator for Diagnostic Tests
#'
#' Calculate the required sample size for testing the hypothesis that sensitivity or specificity of a new test is equal to a particular value.
#'
#' @param theta (numeric): The conjectured sensitivity or specificity.
#' @param theta (numeric): The conjectured sensitivity or specificity.
#' @param alpha (numeric): Significance level (default: 0.05).
#' @param beta (numeric): Type II error rate (default: NULL).
#' @param alternative (character): Type of test, must be one of "less", "greater", or "two.sided" (default: "two.sided").
#' @param R (numeric): Ratio of patients without to with the condition (default: 1).
#' @return An object of class "diag_sample_size_html" containing:
#'   - `sample_size` (list): Required sample size.
#'   - `parameters` (list): Input parameters.
#'   - `html_report` (html): Sample size calculation report.
#' @importFrom stats qnorm
#' @examples
#' n1 <- sample_hypotest_sesp(theta_null=0.5, theta_alter=0.8, alpha=0.05, beta=0.2, alternative="two.sided", R=1) #20
#' print(n1)

#' @export
sample_hypotest_sesp <- function(theta_null, theta_alter, alpha=0.05, beta=NULL, alternative="two.sided", R=1){

  #  (1) ----- validate inputs

  if (!is.numeric(theta_null) || length(theta_null) != 1 || theta_null <= 0 || theta_null >= 1) {
    stop("'theta_null' must be a single numeric value in (0,1).")
  }

  if (!is.numeric(theta_alter) || length(theta_alter) != 1 || theta_alter <= 0 || theta_alter >= 1) {
    stop("'theta_alter' must be a single numeric value in (0,1).")
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value in (0,1).")
  }

  if (!is.null(beta) && (!is.numeric(beta) || length(beta) != 1 || beta <= 0 || beta >= 1)) {
    stop("'beta' must be a single numeric value in (0,1), or NULL.")
  }

  if (!is.character(alternative) || length(alternative) != 1  || !(alternative %in% c("two.sided", "greater", "less"))) {
    stop("'alternative' must be one of 'two.sided', 'greater' or 'less'.")
  }

  if (!is.numeric(R) || length(R) != 1 || R <= 0) {
    stop("'R' must be a single positive numeric value.")
  }

  requireNamespace("htmltools", quietly = TRUE) || stop("Package 'htmltools' is required. Please install it first.")


  #  (2) ----- Calculation modules

  message(paste0("Null hypothesis: θ = ", theta_null))

  if(alternative == "two.sided"){

    message(paste0("Alternative hypothesis: θ ≠ "), theta_null)

  } else if(alternative == "greater") {

    message(paste0("Alternative hypothesis: θ > "), theta_null)

  } else if(alternative == "less") {

    message(paste0("Alternative hypothesis: θ < "), theta_null)
  }

  ## adjust alpha for alternative test
  alpha <- switch(
    alternative,
    "less" = 2 * alpha,
    "greater" = 2 * alpha,
    "two.sided" = alpha,
    stop("Invalid 'alternative' value. Use 'less', 'greater', or 'two.sided'.")
  )


  ## calculate variance function
  V0_theta <- theta_null*(1-theta_null)
  VA_theta <- theta_alter*(1-theta_alter)

  var_function <- c(V0_theta, VA_theta)
  delta <- abs(theta_null - theta_alter)

  ## calculate sample size
  N <- get_diag_sample(
    var_function = var_function,
    alpha = alpha,
    beta = beta,
    delta = delta,
    test_type = "two_diagnostic"
  )

  n_with_condition <- ceiling(N)
  n_without_condition <- ceiling(N*R)
  n_total <- n_with_condition + n_without_condition

  #  (3) ----- Structured outputs

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
          "To test whether the ",
          htmltools::tags$strong(style = paste0("color:", ifelse(theta_alter > 0.5, "#E74C3C", "#27AE60"), ";"),"sensitivity"),
          "or ",
          htmltools::tags$strong(style = paste0("color:", ifelse(theta_alter > 0.5, "#E74C3C", "#27AE60"), ";"),"specificity"),
        " of a new diagnostic test is equal to a pre-determined value"
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
          Parameter = c("Target Accuracy under Null Hypothesis",
                        "Target Accuracy under Alternative Hypothesis",
                        "Confidence Level",
                        "Statistical Power",
                        "Group Allocation"
                      ),
          Symbol    = c("θ0",
                        "θ1",
                        "1 - α",
                        "1 - β",
                        "R"
                      ),
          Value     = c(paste0(round(theta_null*100, 1), "%"),
                        paste0(round(theta_alter*100, 1), "%"),
                        paste0((1-alpha)*100, "%"),
                        paste0((1-beta)*100, "%"),
                        R
                      ),
          Notes     = c("Assumed sensitivity or specificity under null hypothesis",
                        "Assumed sensitivity or specificity under alternative hypothesis",
                        "1 - Type I error rate",
                        "1 - Type II error rate",
                        "Ratio of patients without to with the condition"
                      )
        ),

        options = list(dom = 't', pageLength = 6),
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
            "For testing whether the sensitivity or specificity is equal to a pre-determined value, the null and alternative hypotheses are:"
          ),


          if(alternative == "two.sided") {

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

                "$$ H_{0}: \\theta = \\theta_{0} $$",
                "$$ H_{1}: \\theta ≠ \\theta_{0} $$"
              )
            )

          } else if(alternative == "greater") {

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

                "$$ H_{0}: \\theta = \\theta_{0} $$",
                "$$ H_{1}: \\theta > \\theta_{0} $$"
              )
            )
          } else if(alternative == "less") {

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

                "$$ H_{0}: \\theta = \\theta_{0} $$",
                "$$ H_{1}: \\theta < \\theta_{0} $$"
              )
            )
          },


          htmltools::tags$p(
            "The required sample size is calculated using the formula:"
          ),

          if(alternative == "two.sided") {
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

                "$$ n = \\dfrac{[Z_{\\alpha/2}\\sqrt{V_{0}(\\hat{\\theta})}+Z_{\\beta} \\sqrt{V_{A}(\\hat{\\theta})}]^2}{(\\theta_{0}-\\theta_{1})^2} $$"
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

                "$$ n = \\dfrac{[Z_{\\alpha}\\sqrt{V_{0}(\\hat{\\theta})}+Z_{\\beta} \\sqrt{V_{A}(\\hat{\\theta})}]^2}{(\\theta_{0}-\\theta_{1})^2} $$"
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
            "The variance function in the null hypothsis is calculated through:"
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

              "$$ V_{0}(\\hat{\\theta}) = \\theta_{0} \\times (1-\\theta_{0}) $$"
            )
          ),

          htmltools::tags$p(
            "Similarly, the variance function in the alternative hypothsis is calculated through:"
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

              "$$ V_{A}(\\hat{\\theta}) = \\theta_{1} \\times (1-\\theta_{1}) $$"
            )
          )

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

            if(alternative == "two.sided") {

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
            } else {
              htmltools::tags$div(
                style = "text-align: center; margin: 10px 0;",
                htmltools::tags$p(
                  style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                  sprintf(
                    "$$z_{\\alpha} = z_{%.3f} = %.3f $$",
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

            htmltools::tags$div(
              style = "text-align: center; margin: 10px 0;",
              htmltools::tags$p(
                style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                sprintf(
                  "$$V_{0}(\\hat{\\theta}) = %.3f × (1 - %.3f) = %.3f $$",
                  theta_null, theta_null, theta_null*(1-theta_null)
                )
              )
            ),

            htmltools::tags$div(
              style = "text-align: center; margin: 10px 0;",
              htmltools::tags$p(
                style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                sprintf(
                  "$$V_{A}(\\hat{\\theta}) = %.3f × (1 - %.3f) = %.3f $$",
                  theta_alter, theta_alter, theta_alter*(1-theta_alter)
                )
              )
            ),


            htmltools::tags$p(
              "Now apply the full formula:"
            ),

          htmltools::tags$div(
            style = "text-align: center; margin: 15px 0;",

            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",

              htmltools::HTML(
                sprintf(
                  "$$n = \\dfrac{[(%.3f) \\times \\sqrt{%.3f} + (%.3f) \\times \\sqrt{%.3f}]^2}{(%.3f-%.3f)^2} ≈ %d$$",
                  qnorm(alpha/2),
                  theta_null*(1-theta_null),
                  qnorm(beta),
                  theta_alter*(1-theta_alter),
                  theta_null,
                  theta_alter,
                  N
                )
              )

            )

          ),

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
          paste0(n_with_condition," / ", n_without_condition)
        ),

        htmltools::tags$p(
          paste0("This provides ", (1-alpha)*100, "% confidence with ", (1-beta)*100, "% power that the sensitivity or specificity is equal to a pre-determined value")
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
          "Sample size calculations for hypothesis testing of sensitivity or specificity"
        ),

        ## 2.Assumptions
        htmltools::tags$p(
          style = "margin: 15px 0px; text-indent: -1em; padding-left: 1.5em;",
          htmltools::tags$strong("2. Assumptions:"),
          htmltools::tags$br(),
          "Variance of sensitivity or specificity does not depend on the required sample size"
        ),

        ## 3.Limitations
        htmltools::tags$p(
          style = "margin: 15px 0px; text-indent: -1em; padding-left: 1.5em;",
          htmltools::tags$strong("3. Limitations:"),
          htmltools::tags$br(),
          "Work for sensitivities (or specificities) that are not close to one"
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
          "[3] Arkin CF. How Many Patients Are Necessary to Assess Test Performance? ",
          htmltools::tags$em("JAMA: the Journal of the American Medical Association."),
          " 1990;263(2):275. ",
          htmltools::tags$a(
            href = "doi:https://doi.org/10.1001/jama.1990.03440020109043",
            target = "_blank",
            style = "color: #1a5fb4; text-decoration: none; border-bottom: 1px solid #c1dbf7;
             transition: border-color 0.3s, color 0.3s;",
            onmouseover = "this.style.borderBottomColor='#1a5fb4'; this.style.color='#0d4e9b';",
            onmouseout = "this.style.borderBottomColor='#c1dbf7'; this.style.color='#1a5fb4';",
            "doi:https://doi.org/10.1001/jama.1990.03440020109043"
          ),

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
            "All input parameters (e.g., sensitivity, alpha, beta) should be grounded in real-world data, prior studies, or authoritative literature. We recommend evaluating the impact of ",
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
            "to evaluate the primary diagnostic performance metrics (e.g., sensitivity and specificity) under the specified assumptions."
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
        parameters = list(theta0 = theta_null, theta1 = theta_alter, alpha = alpha, beta = beta, alternative = alternative, R = R),
        html_report = html_report
      ),
      class = "diag_sample_size_html"
    )
  )

}





