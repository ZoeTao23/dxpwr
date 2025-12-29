#' Sample Size Calculator for Diagnostic Tests
#'
#' Calculate the required sample size for comparing sensitivity/specificity between two diagnostic tests.
#'
#' @param theta1 (numeric). The expected sensitivity or specificity of test 1.
#' @param theta2 (numeric). The expected sensitivity or specificity of test 2.
#' @param prob_t1pos_give_t2pos (numeric). Conditional probability of test 1 being positive when test 2 being positive.
#' @param alpha (numeric). Significance level (default: 0.05).
#' @param beta (numeric). Type II error rate.
#' @param paired (logical). Whether the design is paired (`TRUE`) or unpaired (`FALSE`).
#' @return An object of class "diag_sample_size_html" containing:
#'   - `sample_size` (list): Required sample size.
#'   - `parameters` (list): Input parameters.
#'   - `html_report` (html): Sample size calculation report.
#' @examples
#' n <- sample_compare_sesp(theta1=0.8, theta2=0.7, prob_t1pos_give_t2pos=0, alpha=0.05, beta=0.2, paired=FALSE)
#' print(n)
#' @export
sample_compare_sesp <- function(theta1, theta2, prob_t1pos_give_t2pos=NULL, alpha=0.05, beta, paired=FALSE) {


  #  (1) ----- validate inputs

  if (!is.numeric(theta1) || length(theta1) != 1 || theta1 <= 0 || theta1 >= 1) {
    stop("'theta1' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(theta2) || length(theta2) != 1 || theta2 <= 0 || theta2 >= 1) {
    stop("'theta2' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(prob_t1pos_give_t2pos) || length(prob_t1pos_give_t2pos) != 1 || prob_t1pos_give_t2pos < 0 || prob_t1pos_give_t2pos > 1) {
    stop("'prob_t1pos_give_t2pos' must be a single numeric value in [0, 1].")
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(beta) || length(beta) != 1 || beta <= 0 || beta >= 1) {
    stop("'beta' must be a single numeric value in (0, 1).")
  }

  if (!is.logical(paired) || length(paired) != 1) {
    stop("'paired' must be a single logical value TRUE or FALSE.")
  }


  #  (2) -----  calculation modules

  message("Null hypothesis: θ1 = θ2")
  message("Alternative hypothesis: θ1 ≠ θ2")

  ## validate P(T1=1|T2=1)
  if(!paired) {

    prob_t1pos_give_t2pos = theta1
    message("For unpaired design, P(T1=1|T2=1)=P(T1=1)=θ1")

  } else{

    message("For paired design, P(T1=1|T2=1)=1 means a perfect correlation between the test results; P(T1=1|T2=1)=θ1 means a zero correlation between the test results")
  }

  ## calculate variance function
  fi <- theta1 + theta2 - 2 * theta2 * prob_t1pos_give_t2pos
  delta1 <- theta1 - theta2

  V0_delta_theta <- fi
  VA_delta_theta <- fi - delta1^2

  var_function <- c(V0_delta_theta, VA_delta_theta)

  ## calculate sample size
  N <- get_diag_sample(
    var_function = var_function,
    alpha = alpha,
    beta = beta,
    delta = delta1,
    test_type = "two_diagnostic"
  )

  ## calculate number of patients for each test
  if(paired){

    n_test1 <- n_test2 <- ceiling(N)
    n_total <- n_test1

  } else {

    n_test1 <- n_test2 <- ceiling(N)
    n_total <- n_test1 + n_test2
  }


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
        "To compare the ",
        htmltools::tags$strong(style = paste0("color:", ifelse(theta1 > 0.5, "#E74C3C", "#27AE60"), ";"),"sensitivity"), "or",
        htmltools::tags$strong(style = paste0("color:", ifelse(theta1 > 0.5, "#E74C3C", "#27AE60"), ";"),"specificity"),
        "between two diagnostic tests"
      ),

      ## 2.Type
      htmltools::tags$p(
        htmltools::tags$strong("Type:"),
        paste0("Clinical performance evaluation", ifelse(paired, "- paired design", "- unpaired design"))
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
        Parameter = c("Target Accuracy for test 1",
                      "Target Accuracy for test 2",
                      "Conditional Probability",
                      "Confidence Level",
                      "Statistical Power",
                      "Difference in Theta"
        ),
        Symbol    = c("θ1",
                      "θ2",
                      "P(T1=1|T2=1)",
                      "1 - α",
                      "1 - β",
                      "Δ1"
        ),
        Value     = c(theta1,
                      theta2,
                      prob_t1pos_give_t2pos,
                      paste0((1-alpha)*100, "%"),
                      paste0((1-beta)*100, "%"),
                      abs(theta2-theta1)
        ),
        Notes     = c("The expected sensitivity or specificity of test 1 in alternative hypothsis",
                      "The expected sensitivity or specificity of test 2 in alternative hypothsis",
                      "Conditional probability of test 1 being positive when test 2 being positive",
                      "1 - Type I error rate",
                      "1 - Type II error rate",
                      "Difference between sensitivity or specificity under the alternative hypothesis"
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
          "For comparing sensitivity or specificity with a given precision, the null and alternative hypotheses are:"
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

            "$$ H_{0}: \\theta_{1} = \\theta_{2} $$",
            "$$ H_{1}: \\theta_{1} ≠ \\theta_{2} $$"
          )
        ),

        htmltools::tags$p(
          "The required sample size is calculated using the formula:"
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

            "$$ n = \\dfrac{[Z_{\\alpha/2} \\sqrt{V_{0}(\\hat{\\theta}_{1}-\\hat{\\theta}_{2})} + Z_{\\beta} \\sqrt{V_{A}(\\hat{\\theta}_{1}-\\hat{\\theta}_{2})}]^2}{(Δ_{1})^2} $$"
          )
        )

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

            paste0("$$ V_{0}(\\hat{\\theta}_{1}-\\hat{\\theta}_{2}) = \\theta_{1} + \\theta_{2} - 2 \\times \\theta_{2} \\times P(T_{1}=1|T_{2}=1) $$")
          )
        ),


        htmltools::tags$p(
          "The variance function in the alternative hypothsis is calculated through:"
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

            paste0("$$ V_{A}(\\hat{\\theta}_{1}-\\hat{\\theta}_{2}) = V_{0}(\\hat{\\theta}_{1}-\\hat{\\theta}_{2}) - Δ_{1}^2 $$")
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
        ),


        htmltools::tags$p(
          "Next, compute the variance function term:"
        ),

        htmltools::tags$div(
          style = "text-align: center; margin: 10px 0;",
          htmltools::tags$p(
            style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
            sprintf(
              paste0("$$ V_{0}(\\hat{\\theta}_{1}-\\hat{\\theta}_{2}) = %.3f + %.3f - 2 \\times %.3f \\times %.3f = %.3f $$"),
              theta1, theta2, theta2, prob_t1pos_give_t2pos, V0_delta_theta
            )
          )
        ),

        htmltools::tags$div(
          style = "text-align: center; margin: 10px 0;",
          htmltools::tags$p(
            style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
            sprintf(
              paste0("$$ V_{A}(\\hat{\\theta}_{1}-\\hat{\\theta}_{2}) = %.3f - (%.3f - %.3f)^2= %.3f $$"),
              V0_delta_theta, theta1, theta2, VA_delta_theta
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
                "$$n = \\dfrac{[(%.3f)\\times \\sqrt{%.3f}+(%.3f)\\times \\sqrt{%.3f}]^2}{(%.3f-%.3f)^2} ≈ %d$$",
                qnorm(alpha/2),
                V0_delta_theta,
                qnorm(beta),
                VA_delta_theta,
                theta1,
                theta2,
                N
              )
            )

          )

        ),


        htmltools::tags$p(
          "Thus,"
        ),


        if(paired) {
          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",
            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
              sprintf(
                paste0("$$ N_{1} = N_{2} = n = %d $$"),
                n_total
              )
            )
          )

        } else {
          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",
            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
              sprintf(
                paste0("$$ N_{1} = n = %d, N_{2} = n = %d$$"),
                n_test1, n_test2
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
        ceiling(n_total)
      ),

      htmltools::tags$h4(style = "margin-top: 0;", "Participants for test 1 / test 2"),

      htmltools::tags$p(
        style = "font-size: 24px; font-weight: bold; color: #2E86C1;",
        paste0(n_test1," / ", n_test2)
      ),

      htmltools::tags$p(
        paste0("This provides ", (1-beta)*100, "% power"),
        paste0(" with ", (1-alpha)*100, "% confidence"),
        "to detect a difference in sensitivity or specificity"
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
        "Sample size calculations for comparing sensitivity or specificity"
      ),

      ## 2.Assumptions


        htmltools::tags$p(
          style = "margin: 15px 0px; text-indent: -1em; padding-left: 1.5em;",
          htmltools::tags$strong("2. Assumptions:"),
          htmltools::tags$br(),
          "The test results between patients must be mutually independent"

        ),


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
        "[1] Connor RJ. Sample Size for Testing Differences in Proportions for the Paired-Sample Design. ",
        htmltools::tags$em("Biometrics."),
        " 1987;43(1):207. ",
        htmltools::tags$a(
          href = "doi:https://doi.org/10.2307/2531961",
          target = "_blank",
          style = "color: #1a5fb4; text-decoration: none; border-bottom: 1px solid #c1dbf7;
             transition: border-color 0.3s, color 0.3s;",
          onmouseover = "this.style.borderBottomColor='#1a5fb4'; this.style.color='#0d4e9b';",
          onmouseout = "this.style.borderBottomColor='#c1dbf7'; this.style.color='#1a5fb4';",
          "doi:https://doi.org/10.2307/2531961"
        ),


        htmltools::tags$br(),
        "[2] Beam CA. Strategies for Improving Power in Diagnostic Radiology research. ",
        htmltools::tags$em("American Journal of Roentgenology."),
        " 1992;159(3):631-637. ",
        htmltools::tags$a(
          href = "doi:https://doi.org/10.2214/ajr.159.3.1503041",
          target = "_blank",
          style = "color: #1a5fb4; text-decoration: none; border-bottom: 1px solid #c1dbf7;
             transition: border-color 0.3s, color 0.3s;",
          onmouseover = "this.style.borderBottomColor='#1a5fb4'; this.style.color='#0d4e9b';",
          onmouseout = "this.style.borderBottomColor='#c1dbf7'; this.style.color='#1a5fb4';",
          "doi:https://doi.org/10.2214/ajr.159.3.1503041"
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
          "All input parameters (e.g., theta, alpha, beta) should be grounded in real-world data, prior studies, or authoritative literature. We recommend evaluating the impact of ",
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
          "to compare the primary diagnostic performance metrics (e.g., sensitivity and specificity) in ",
          ifelse(paired, "a paired study design " , "an unpaired study design "),
          "under the specified assumptions."
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
        sample_size = list(n_total = n_total, n_test1 = n_test1, n_test2 = n_test2),
        parameters = list(theta1 = theta1, theta2 = theta2, prob_t1pos_give_t2pos = prob_t1pos_give_t2pos, alpha = alpha, beta = beta, paired = paired),
        html_report = html_report
      ),
      class = "diag_sample_size_html"
    )
  )


}
