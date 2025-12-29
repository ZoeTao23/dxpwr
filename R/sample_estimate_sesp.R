#' Sample Size Calculator for Diagnostic Tests
#'
#' Calculate the required sample size for estimating sensitivity or specificity of a single diagnostic test.
#'
#' @param theta (numeric). The conjectured sensitivity or specificity.
#' @param alpha (numeric). Significance level (default: 0.05).
#' @param beta (numeric). Type II error rate (default: NULL).
#' @param L (numeric). The desired length of one-half of the (1-α)×100% confidence interval for sensitivity/specificity.
#' @return An object of class "diag_sample_size_html" containing:
#'   - `sample_size`: Required sample size.
#'   - `parameters`: Input parameters.
#'   - `html_report`: Sample size calculation report.
#' @importFrom stats qnorm
#' @examples
#' n1 <- sample_estimate_sesp(theta=0.8, alpha=0.05, beta=0.2, L=0.05) #502
#' print(n1)
#' n2 <- sample_estimate_sesp(theta=0.8, alpha=0.05, L=0.05)
#' print(n2)
#' @export
sample_estimate_sesp <- function(theta, alpha=0.05, beta=NULL, L){

  #  (1) ----- validate inputs

  if (!is.numeric(theta) || length(theta) != 1 || theta <= 0 || theta >= 1) {
    stop("'theta' must be a single numeric value in (0,1).")
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value in (0,1).")
  }

  if (!is.null(beta) && (!is.numeric(beta) || length(beta) != 1 || beta <= 0 || beta >= 1)) {
    stop("'beta' must be a single numeric value in (0,1), or NULL.")
  }

  if (!is.numeric(L) || length(L) != 1 || L <= 0 || L > 0.5 ) {
    stop("'L' must be a single numeric value in (0, 0.5).")
  }

  requireNamespace("htmltools", quietly = TRUE) || stop("Package 'htmltools' is required. Please install it first.")


  #  (2) ----- Calculation modules

  ## calculate variance function
  var_function <- theta*(1-theta)

  ## calculate sample size
  N <- get_diag_sample(
    var_function = var_function,
    alpha = alpha,
    beta = beta,
    delta = L,
    test_type = "one_diagnostic"
  )

  n_total <- N

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
          "To evaluate the ",
          htmltools::tags$strong(style = paste0("color:", ifelse(theta > 0.5, "#E74C3C", "#27AE60"), ";"),"sensitivity"),
          "or ",
          htmltools::tags$strong(style = paste0("color:", ifelse(theta > 0.5, "#E74C3C", "#27AE60"), ";"),"specificity"),
        " of a new diagnostic test"
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
                        "Confidence Level",
                        "Statistical Power",
                        "Margin of Error"
                      ),
          Symbol    = c("θ",
                        "1 - α",
                        "1 - β",
                        "± L"
                      ),
          Value     = c(paste0(round(theta*100, 1), "%"),
                        paste0((1-alpha)*100, "%"),
                        ifelse(is.null(beta), "Not specified", paste0((1-beta)*100, "%")),
                        paste0("±", L*100, "%")
                      ),
          Notes     = c("Desired sensitivity or specificity",
                        "1 - Type I error rate",
                        "1 - Type II error rate",
                        "Half-width of CI"
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
              "For estimating sensitivity/specificity with a given precision, the required sample size is calculated using the formula:"
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

                  "$$ n = \\dfrac{[Z_{\\alpha/2} \\sqrt{V(\\hat{\\theta})}]^2}{L^2} $$"
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

                  "$$ n = \\dfrac{[(Z_{\\alpha/2}+Z_{\\beta}) \\sqrt{V(\\hat{\\theta})}]^2}{L^2} $$"
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

                "$$ V(\\hat{\\theta}) = \\theta \\times (1-\\theta) $$"
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

              htmltools::tags$div(
                style = "text-align: center; margin: 10px 0;",
                htmltools::tags$p(
                  style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                  sprintf(
                    "$$V(\\hat{\\theta}) = %.3f × (1 - %.3f) = %.3f $$",
                    theta, theta, theta*(1-theta)
                  )
                )
              ),


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
                        theta*(1-theta),
                        L,
                        n_total
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
                        "$$n = \\dfrac{\\{[(%.3f)+(%.3f)]\\times \\sqrt{%.3f}\\}^2}{(%.3f)^2} ≈ %d$$",
                        qnorm(alpha/2),
                        qnorm(beta),
                        theta*(1-theta),
                        L,
                        n_total
                      )
                    )

                  )

                )
              }

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
          n_total
        ),

        htmltools::tags$p(
          paste0("This provides ", (1-alpha)*100, "% confidence that the estimated sensitivity or specificity"),
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
          "Sample size calculations for estimating sensitivity or specificity"
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
        sample_size = list(n_total = n_total),
        parameters = list(theta = theta, alpha = alpha, beta = beta, L = L),
        html_report = html_report
      ),
      class = "diag_sample_size_html"
    )
  )

}





