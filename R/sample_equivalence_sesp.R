#' Sample Size Calculator for Diagnostic Tests
#'
#' Calculate the required sample size for testing whether the sensitivity/specificity of a new test is equivalent to an existing test.
#'
#' @param theta_S (numeric): The expected sensitivity (or specificity) of an existing test (standard test).
#' @param theta_E (numeric): The expected sensitivity (or specificity) of a new test (experimental test).
#' @param prob_t1pos_give_t2pos (numeric): Conditional probability of the existing test being positive when the new test being positive.
#' @param alpha (numeric): Significance level (default: 0.05).
#' @param beta (numeric): Type II error rate.
#' @param delta (numeric): The prespecified upper clinical limits for equivalence.
#' @param paired (character): Whether the design is paired (`TRUE`) or unpaired (`FALSE`).
#' @return An object of class "diag_sample_size_html" containing:
#'   - `sample_size` (list): Required sample size.
#'   - `parameters` (list): Input parameters.
#'   - `html_report` (html): Sample size calculation report.
#' @examples
#' n <- sample_equivalence_sesp(theta_S=0.9, theta_E=0.88, prob_t1pos_give_t2pos=0, delta=0.05, alpha=0.05, beta=0.2, R=1, paired=FALSE)
#' print(n)
#' @export
sample_equivalence_sesp <- function(theta_S, theta_E, prob_t1pos_give_t2pos=NULL, delta, alpha=0.05, beta, R=1, paired=TRUE) {


  #  (1) ----- validate inputs

  if (!is.numeric(theta_S) || length(theta_S) != 1 || theta_S <= 0 || theta_S >= 1) {
    stop("'theta_S' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(theta_E) || length(theta_E) != 1 || theta_E <= 0 || theta_E >= 1) {
    stop("'theta_E' must be a single numeric value in (0, 1).")
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

  if (!is.numeric(delta) || length(delta) != 1 || delta < 0 || delta > 1) {
    stop("'delta' must be a single numeric value in [0, 1].")
  }

  if (!is.logical(paired) || length(paired) != 1) {
    stop("'paired' must be a single logical value TRUE or FALSE.")
  }


  #  (2) -----  calculation modules

  message(paste0("Null hypothesis: θs - θe ≤ ", -delta, " or θs - θe ≥ ", delta))
  message(paste0("Alternative hypothesis: ", -delta, " < θs - θe < ", delta))

  ## validate P(T1=1|T2=1)
  if(!paired) {

    prob_t1pos_give_t2pos = theta_S
    message("For unpaired design, P(Ts=1|Te=1)=P(Ts=1)=θs")

  } else{

    message("For paired design, P(Ts=1|Te=1)=1 means a perfect correlation between the test results; P(Ts=1|Te=1)=θs means a zero correlation between the test results")
  }



  ## calculate variance function
  V0_delta_theta <- theta_S + theta_E - 2 * theta_E * prob_t1pos_give_t2pos
  VA_delta_theta <- V0_delta_theta - (theta_S - theta_E)^2

  var_function <- VA_delta_theta


  ## calculate delta
  if(theta_S - theta_E > 0) {

    delta1 <- abs(delta - (theta_S - theta_E))

  } else if (theta_S - theta_E < 0) {

    delta1 <- abs(delta + (theta_S - theta_E))

  } else if (theta_S - theta_E == 0){

    beta = beta/2
    delta1 <- delta

  }

  ## calculate sample size
  N <- get_diag_sample(
    var_function = var_function,
    alpha = alpha,
    beta = beta,
    delta = delta1,
    test_type = "equivalence"
  )


  ## calculate number of patients for each test
  if (paired) {

    n_test1_with_condition <- ceiling(N)
    n_test1_without_condition <- ceiling(R*N)
    n_test1 <- n_test1_with_condition + n_test1_without_condition

    n_test2_with_condition <- ceiling(N)
    n_test2_without_condition <- ceiling(R*N)
    n_test2 <- n_test2_with_condition + n_test2_without_condition

    n_total <- n_test1

  } else {

    n_test1_with_condition <- ceiling(N)
    n_test1_without_condition <- ceiling(R*N)
    n_test1 <- n_test1_with_condition + n_test1_without_condition

    n_test2_with_condition <- ceiling(N)
    n_test2_without_condition <- ceiling(R*N)
    n_test2 <- n_test2_with_condition + n_test2_without_condition

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
        "To test whether the ",
        htmltools::tags$strong(style = paste0("color:", ifelse(theta_S > 0.5, "#E74C3C", "#27AE60"), ";"),"sensitivity"), "or",
        htmltools::tags$strong(style = paste0("color:", ifelse(theta_S > 0.5, "#E74C3C", "#27AE60"), ";"),"specificity"),
        " of a new test is equivalent to an existing test"
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
        Parameter = c("Target Accuracy for the Existing Test",
                      "Target Accuracy for the New Test",
                      "Conditional Probability",
                      "Confidence Level",
                      "Statistical Power",
                      "Equivalence Margin"
        ),
        Symbol    = c("θs",
                      "θe",
                      "P(Ts=1|Te=1)",
                      "1 - α",
                      "1 - β",
                      "(ΔL, ΔU)"
        ),
        Value     = c(theta_S,
                      theta_E,
                      prob_t1pos_give_t2pos,
                      paste0((1-alpha)*100, "%"),
                      paste0((1-beta)*100, "%"),
                      paste0("(",-delta,", ",delta,")")
        ),
        Notes     = c("The expected sensitivity or specificity of the existing test",
                      "The expected sensitivity or specificity of the new test",
                      "Conditional probability of the existing test being positive when the new test being positive",
                      "1 - Type I error rate",
                      "1 - Type II error rate",
                      "The prespecified lower and upper clinical limits for equivalence"
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
          "For the sensitivity or specificity equivalence testing, the null and alternative hypotheses are formulated as:"
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

            "$$ H_{0}: \\theta_{S} - \\theta_{E} \\leq \\Delta_{L} \\quad\\text{or}\\quad \\theta_{S} - \\theta_{E} \\geq \\Delta_{U} $$",
            "$$ H_{1}: \\Delta_{L} < \\theta_{S} - \\theta_{E} < \\Delta_{U} $$"
          )
        ),

        htmltools::tags$p(
          "The required sample size is calculated using the formula:"
        ),

        if (theta_S - theta_E == 0) {

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

              "$$ n = \\dfrac{(Z_{\\alpha}+Z_{\\beta/2})^2 V_{A}(\\hat{\\theta}_{S}-\\hat{\\theta}_{E})}{(Δ_{U})^2} $$"
            )
          )

        } else if(theta_S - theta_E != 0) {

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

              "$$ n = \\dfrac{(Z_{\\alpha}+Z_{\\beta})^2 V_{A}(\\hat{\\theta}_{S}-\\hat{\\theta}_{E})}{(Δ_{U}-|\\theta_{S}-\\theta_{E}|)^2} $$"
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

            paste0("$$ V_{0}(\\hat{\\theta}_{S}-\\hat{\\theta}_{E}) = \\theta_{S} + \\theta_{E} - 2 \\times \\theta_{E} \\times P(T_{S}=1|T_{E}=1) $$")
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

            paste0("$$ V_{A}(\\hat{\\theta}_{S}-\\hat{\\theta}_{E}) = V_{0}(\\hat{\\theta}_{S}-\\hat{\\theta}_{E}) - (\\theta_{S} - \\theta_{E})^2 $$")
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

        if (theta_S - theta_E == 0) {

          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",
            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
              sprintf(
                "$$z_{\\alpha} = z_{%.2f} = %.3f $$",
                alpha, qnorm(alpha)
              ),
              sprintf(
                "$$z_{\\beta/2} = z_{%.2f} = %.3f $$",
                beta, qnorm(beta)
              )
            )
          )

        } else if (theta_S - theta_E != 0) {

          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",
            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
              sprintf(
                "$$z_{\\alpha} = z_{%.3f} = %.3f $$",
                alpha, qnorm(alpha)
              ),
              sprintf(
                "$$z_{\\beta} = z_{%.3f} = %.3f $$",
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
              paste0("$$ V_{0}(\\hat{\\theta}_{S}-\\hat{\\theta}_{E}) = %.3f + %.3f - 2 \\times %.3f \\times %.3f = %.3f $$"),
              theta_S, theta_E, theta_E, prob_t1pos_give_t2pos, V0_delta_theta
            )
          )
        ),

        htmltools::tags$div(
          style = "text-align: center; margin: 10px 0;",
          htmltools::tags$p(
            style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
            sprintf(
              paste0("$$ V_{A}(\\hat{\\theta}_{S}-\\hat{\\theta}_{E}) = %.3f - (%.3f - %.3f)^2= %.3f $$"),
              V0_delta_theta, theta_S, theta_E, VA_delta_theta
            )
          )
        ),


        htmltools::tags$p(
          "Now apply the full formula:"
        ),

        if(theta_S - theta_E == 0) {

          htmltools::tags$div(
            style = "text-align: center; margin: 15px 0;",

            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",

              htmltools::HTML(
                sprintf(
                  "$$n = \\dfrac{[(%.3f) + (%.3f)]^2 \\times %.3f}{(%.3f)^2} ≈ %d$$",
                  qnorm(alpha),
                  qnorm(beta),
                  VA_delta_theta,
                  delta,
                  N
                )
              )

            )

          )
        } else if (theta_S - theta_E != 0) {

          htmltools::tags$div(
            style = "text-align: center; margin: 15px 0;",

            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",

              htmltools::HTML(
                sprintf(
                  "$$n = \\dfrac{[(%.3f) + (%.3f)]^2 \\times %.3f}{(%.3f - |%.3f-%.3f|)^2} ≈ %d$$",
                  qnorm(alpha),
                  qnorm(beta),
                  VA_delta_theta,
                  delta,
                  theta_S,
                  theta_E,
                  N
                )
              )

            )

          )

        },


        htmltools::tags$p(
          "Thus,"
        ),


        if(paired) {

          htmltools::tagList(
            htmltools::tags$div(
              style = "text-align: center; margin: 10px 0;",
              htmltools::tags$p(
                style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                sprintf(
                  paste0("$$ N_{S+} = N_{E+} = n = %d $$"),
                  n_test1_with_condition
                )
              )
            ),

            htmltools::tags$div(
              style = "text-align: center; margin: 10px 0;",
              htmltools::tags$p(
                style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                sprintf(
                  paste0("$$ N_{S-} = N_{E-} = R \\times n = %d $$"),
                  n_test1_without_condition
                )
              )
            )
          )

        } else {

          htmltools::tagList(
            htmltools::tags$div(
              style = "text-align: center; margin: 10px 0;",
              htmltools::tags$p(
                style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                sprintf(
                  paste0("$$ N_{S+} = n = %d, N_{E+} = n = %d$$"),
                  n_test1_with_condition, n_test2_with_condition
                )
              )
            ),

            htmltools::tags$div(
              style = "text-align: center; margin: 10px 0;",
              htmltools::tags$p(
                style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                sprintf(
                  paste0("$$ N_{S-} = R \\times n = %d, N_{E-} = R \\times n = %d$$"),
                  n_test1_without_condition, n_test2_without_condition
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
        "to test equivalence in sensitivity or specificity"
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
        "Sample size calculations for the sensitivity or specificity equivalence testing"
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
        ),


        htmltools::tags$br(),
        "[3] Dongsheng T. Two one-sided Tests Procedures in Establishing Therapeutic Equivalence with Binary Clinical endpoints: Fixed Sample Performances and Sample Size Derformination. ",
        htmltools::tags$em("Journal of Statistical Computation and Simulation."),
        " 1997;59(3):271-290. ",
        htmltools::tags$a(
          href = "doi:https://doi.org/10.1080/00949659708811860",
          target = "_blank",
          style = "color: #1a5fb4; text-decoration: none; border-bottom: 1px solid #c1dbf7;
             transition: border-color 0.3s, color 0.3s;",
          onmouseover = "this.style.borderBottomColor='#1a5fb4'; this.style.color='#0d4e9b';",
          onmouseout = "this.style.borderBottomColor='#c1dbf7'; this.style.color='#1a5fb4';",
          "doi:https://doi.org/10.1080/00949659708811860"
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
          "to demonstrate that the difference in primary diagnostic performance metrics (e.g., sensitivity and specificity) falls within the equivalence margin in ",
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
        sample_size = list(n_total = n_total, n_test1 = n_test1, n_test2 = n_test2, n_test1_with_condition = n_test1_with_condition, n_test1_without_condition = n_test1_without_condition, n_test2_with_condition = n_test2_with_condition, n_test2_without_condition = n_test2_without_condition),
        parameters = list(theta_S = theta_S, theta_E = theta_E, prob_t1pos_give_t2pos = prob_t1pos_give_t2pos, alpha = alpha, beta = beta, paired = paired),
        html_report = html_report
      ),
      class = "diag_sample_size_html"
    )
  )

}
