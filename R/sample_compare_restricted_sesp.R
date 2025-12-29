#' Sample Size Calculator for Diagnostic Tests
#'
#' Calculate the required sample size for comparing sensitivity at a fixed false positive rate (FPR) or specificity at a fixed false negative rate (FNR) of two diagnostic tests.
#'
#' @param metric (character): Performance metric to be evaluated, must be either "sens" or "sepc":
#'   - `"sens"`: Sensitivity at the fixed FPR.
#'   - `"spec"`: Specificity at the fixed FNR.
#' @param A_null (numeric): AUC under the null hypothesis.
#' @param b_null (numeric): Ratio of the standard deviations of the distributions of test results for patients without versus with the condition, σ_without/σ_with (default: 1).
#' @param A1_alter (numeric): AUC of test 1 under the alternative hypothesis.
#' @param b1_alter (numeric): Ratio of standard deviations of test 1 under the alternative hypothesis.
#' @param A2_alter (numeric): AUC of test 2 under the alternative hypothesis.
#' @param b2_alter (numeric): Ratio of standard deviations of test 2 under the alternative hypothesis.
#' @param e (numeric): For "sens" provide false positive rate value, and for "spec" provide false negative rate value.
#' @param rD (numeric): The correlation of the underlying bivariate binormal distribution for patients with the condition (default: 0.5).
#' @param rN (numeric): The correlation of the underlying bivariate binormal distribution for patients without the condition (default: 0.5).
#' @param alpha (numeric): Significance level (default: 0.05).
#' @param beta (numeric): Type II error rate.
#' @param paired (logical). Whether the design is paired (`TRUE`) or unpaired (`FALSE`).
#' @param dist (Character): The assumption for choosing the variance function (default: "binorm"):
#'   - `"binorm"`: Assume that the unobserved, underlying test results follow a binormal distribution, but the observed test results, either continuous or ordinal, do not necessarily have a binormal distribution.
#'   - `"obs_binorm"`: Assume that the observed test results are on a truly continuous scale and they follow a binormal distribution (or can be transformed to a binormal distribution).
#' @param R (numeric): Ratio of patients with and without the condition (default: 1).
#' @return An object of class "diag_sample_size_html" containing:
#'   - `sample_size` (list): Required sample size.
#'   - `parameters` (list): Input parameters.
#'   - `html_report` (html): Sample size calculation report.
#' @examples
#' n <- sample_compare_restricted_sesp(metric="sens", A_null=0.8, b_null=1, A1_alter=0.8, b1_alter=1, A2_alter=0.9, b2_alter=1, e=0.1, rN=0.5, rD=0.5, alpha=0.05, beta=0.2, paired=TRUE, R=1, dist ="binorm")
#' print(n)
#' @export
sample_compare_restricted_sesp <- function(metric, A_null, b_null=1, A1_alter, b1_alter=1, A2_alter, b2_alter=1, e, rN=0.5, rD=0.5, alpha=0.05, beta, paired=TRUE, R=1) {

  #  (1) ----- validate inputs
  if (length(metric) != 1 || !metric %in% c("sens", "spec")) {
    stop("'metric' must be either 'sens' or 'spec'.")
  }

  if (!is.numeric(A_null) || length(A_null) != 1 || A_null <= 0 ||A_null >= 1) {
    stop("'A_null' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(A1_alter) || length(A1_alter) != 1 || A1_alter <= 0 || A1_alter >= 1) {
    stop("'A1_alter' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(A2_alter) || length(A2_alter) != 1 || A2_alter <= 0 || A2_alter >= 1) {
    stop("'A2_alter' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(b_null) || length(b_null) != 1 || b_null <= 0) {
    stop("'b_null' must be a single positive numeric value.")
  }

  if (!is.numeric(b1_alter) || length(b1_alter) != 1 || b1_alter <= 0) {
    stop("'b1_alter' must be a single positive numeric value.")
  }

  if (!is.numeric(b2_alter) || length(b2_alter) != 1 || b2_alter <= 0) {
    stop("'b2_alter' must be a single positive numeric value.")
  }

  if (!is.numeric(e) || length(e) != 1 || e <= 0 || e >= 1) {
    stop("'e' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(rD) || length(rD) != 1 || rD < -1 || rD > 1) {
    stop("'rD' must be a single numeric value in [-1, 1].")
  }

  if (!is.numeric(rN) || length(rN) != 1 || rN < -1 || rN > 1) {
    stop("'rN' must be a single numeric value in [-1, 1].")
  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(beta) || length(beta) != 1|| beta <= 0 || beta >= 1) {
    stop("beta must be a single numeric value in (0, 1).")
  }

  if (!is.logical(paired) || length(paired) != 1) {
    stop("'paired' must be a single logical value TRUE or FALSE.")
  }

  if (!is.numeric(R) || length(R) != 1 || R <= 0) {
    stop("'R' must be a single positive numeric value.")
  }

  if (!is.character(dist) || length(dist) != 1 || !(dist %in% c("binorm", "obs_binorm"))) {
    stop("'dist' must be either 'binorm' or 'obs_binorm'.")
  }

  #  (2) ----- calculation modules

  ## assign empirical correlation
  if(paired){

    message("For paired design, rD = rN = 0.5 is a typical value for the correlation between tests")
    message("If not specified, rD and rN default to 0.5")

  } else {

    message("For unpaired design, rD = rN = 0")
    r <- 0
    rN <- rD <- r
  }


  ## calcuate variance function
  a_null <- get_binorm_params(A_null, b_null)$a
  a1_alter <- get_binorm_params(A1_alter, b1_alter)$a
  a2_alter <- get_binorm_params(A2_alter, b2_alter)$a

  z1_alter <- a1_alter + b1_alter * qnorm(e)
  z2_alter <- a2_alter + b2_alter * qnorm(e)

  if(dist == "binorm") {

    V_theta <- function(a, b, e, R) {

      g <- qnorm(e)
      V <- 1 + b^2 / R + a^2 / 2 + g^2 * b^2 * (1 + R) / (2 * R)

      return(V)
    }

  } else if(dist == "obs_binorm") {

    V_theta <- function(a, b, e, R) {

      g <- qnorm(e)
      V <- 1 + b^2 / R + a^2 / 2 + g^2 * b^2 * (1 + R) / (2 * R) + g * a * b

      return(V)
    }

  } else {

    stop("Invalid 'dist' value. Use 'binorm' or 'obs_binorm'.")

  }


  C_delta_theta <- function(a1, a2, b1, b2, e, rD, rN) {

    g <- qnorm(e)
    C <- rD + (rN * b1 * b2) / R + (rD^2 * a1 * a2) / 2 + g ^ 2 * b1 * b2 * (rN^2 + R * rD^2) / (2 * R) + g* rD^2 * (a1 * b2 + a2 * b1) / 2

    return(C)
  }

  V0_A1 <- V0_A2 <- V_theta(a_null, b_null, e, R)
  VA_A1 <- 1 + V_theta(a1_alter, b1_alter, e, R)
  VA_A2 <- 1 + V_theta(a2_alter, b2_alter, e, R)

  C0_delta_A <- C_delta_theta(a_null, a_null, b_null, b_null, e, rD, rN)
  CA_delta_A <- C_delta_theta(a1_alter, a2_alter, b1_alter, b2_alter, e, rD, rN)

  V0_delta_A <- V0_A1 + V0_A2 - 2*C0_delta_A
  VA_delta_A <- VA_A1 + VA_A2 - 2*CA_delta_A
  var_function <- c(V0_delta_A, VA_delta_A)

  ## calculate delta
  delta1 <- abs(z1_alter - z2_alter)


  ## calculate sample size with and without condition
  N <- get_diag_sample(
    var_function = var_function,
    alpha = alpha,
    beta = beta,
    delta = delta1,
    test_type = "two_diagnostic"
  )

  if(paired){

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
        "To compare the ",
        htmltools::tags$strong(style = paste0("color:", ifelse(a1_alter > 0.5, "#E74C3C", "#27AE60"), ";"), ifelse(metric=="sens", "sensitivity at fixed FPR", "specificity at fixed FNR")),
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
        Parameter = c("Target Accuracy in Null Hypothesis",
                      "Target Accuracy for Test 1",
                      "Target Accuracy for Test 2",
                      "Binominal Parameter in Null Hypothesis",
                      "Binominal Parameter for Test 1",
                      "Binominal Parameter for Test 2",
                      ifelse(metric=="sens", "False Positive Rate", "False Negative Rate"),
                      "Correlation for Patients with Condition",
                      "Correlation for Patients without Condition",
                      "Confidence Level",
                      "Statistical Power",
                      "Group Allocation"

        ),
        Symbol    = c("A_null",
                      "A1_alter",
                      "A2_alter",
                      "(a_null, b_null)",
                      "(a1_alter, b1_alter)",
                      "(a2_alter, b2_alter)",
                      "e",
                      "rN",
                      "rD",
                      "1 - α",
                      "1 - β",
                      "R"
        ),
        Value     = c(ifelse(is.null(A_null),"Not specified", A_null),
                      A1_alter,
                      A2_alter,
                      ifelse(is.null(A_null),"Not specified", paste0("(", round(a_null,3), ", ", round(b_null,3), ")")),
                      paste0("(", round(a1_alter,3), ", ", round(b1_alter,3), ")"),
                      paste0("(", round(a2_alter,3), ", ", round(b2_alter,3), ")"),
                      e,
                      ifelse(is.na(rN), "Not specified", rN),
                      ifelse(is.na(rD), "Not specified", rD),
                      paste0((1-alpha)*100, "%"),
                      paste0((1-beta)*100, "%"),
                      R
        ),
        Notes     = c("The expected AUC in null hypothesis",
                      "The expected AUC of test 1 in alternative hypothesis",
                      "The expected AUC of test 2 in alternative hypothesis",
                      "a = (μ0_with - μ0_without)/σ0_with; b = σ0_without/σ0_with, to determine the shape of ROC in null hypothesis",
                      "a = (μ1_with - μ1_without)/σ1_with; b = σ1_without/σ1_with, to determine the shape of ROC of test 1 in alternative hypothesis",
                      "a = (μ2_with - μ2_without)/σ2_with; b = σ2_without/σ2_with, to determine the shape of ROC of test 2 in alternative hypothesis",
                      paste0("Constraint on ", ifelse(metric=="sens", "sensitivity", "specificity")),
                      "The correlation of the underlying bivariate binormal distribution for patients with the condition",
                      "The correlation of the underlying bivariate binormal distribution for patients without the condition",
                      "1 - Type I error rate",
                      "1 - Type II error rate",
                      "Ratio of patients with and without the condition"
        )
      ),

      options = list(dom = 't', pageLength = 13),
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

        if(metric=="sens") {

          htmltools::tagList(
            htmltools::tags$p(
              "For comparing sensitivity with a given precision, the null and alternative hypotheses are:"
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

                "$$ H_{0}: Se_{1,FPR=e} = Se_{2,FPR=e} $$",
                "$$ H_{1}: Se_{1,FPR=e} ≠ Se_{2,FPR=e} $$"
              )
            )
          )
        }
        else{

          htmltools::tagList(
            htmltools::tags$p(
              "For comparing specificity with a given precision, the null and alternative hypotheses are:"
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

                "$$ H_{0}: Sp_{1,FNR=e} = Sp_{2,FNR=e} $$",
                "$$ H_{1}: Sp_{1,FNR=e} ≠ Sp_{2,FNR=e} $$"
              )
            )
          )
        },


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

            "$$ n = \\dfrac{[Z_{\\alpha/2} \\sqrt{V_{0}(\\hat{A}_{1}-\\hat{A}_{2})} + Z_{\\beta} \\sqrt{V_{A}(\\hat{A}_{1}-\\hat{A}_{2})}]^2}{(Δ_{1})^2} $$"
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

        htmltools::tagList(

          htmltools::tags$p(
            "The variance function of the estimated difference in AUC is calculated through:"
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

              "$$ V(\\hat{A}_{1}-\\hat{A}_{2}) = V(\\hat{A}_{1}) + V(\\hat{A}_{2}) - 2C(\\hat{A}_{1},\\hat{A}_{2}) $$"
            )
          ),

          htmltools::tags$p(
            "Where, the variance fuction of AUC is"
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

                paste0("$$ \\hat{V}(z[\\hat{Se}_{",ifelse(metric=="sens","FPR","FNR"),"=e}]) = 1+b^2/R+a^2/2+b^2[\\Phi^{-1}(e)]^2(1+R)/(2R) $$")
              )
            )

          } else if (dist=="obs_binorm") {

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

                paste0("$$ \\hat{V}(z[\\hat{Se}_{",ifelse(metric=="sens","FPR","FNR"),"=e}]) = 1+b^2/R+a^2/2+b^2[\\Phi^{-1}(e)]^2(1+R)/(2R)+ab\\Phi^{-1}(e) $$")
              )
            )

          },

          htmltools::tags$p(
            paste0(ifelse(paired, "And for paired design, ", "And for unpaired design, "), "the covariance function is:")
          ),

          if(paired) {

            htmltools::tagList(
              htmltools::tags$div(
                style = "text-align: center; margin: 10px 0;",
                htmltools::tags$p(
                  style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                  paste0(ifelse(metric=="sens","$$\\hat{C}((\\hat{Se}_{FPR=e})_{1}, (\\hat{Se}_{FPR=e})_{2})","$$\\hat{C}((\\hat{Sp}_{FNR=e})_{1}, (\\hat{Sp}_{FNR=e})_{2})"),"= (r_{D}+\\frac{r_{N}b_{1}b_{2}}{R}+\\frac{r_{D}^2a_{1}a_{2}}{2}) + \\frac{g^2b_{1}b_{2}(r_{N}^2+Rr_{D}^2)}{2R} + \\frac{gr_{D}^2(a_{1}b_{2}+a_{2}b_{1})}{2} \\times ()  $$")

                )
              )

            )

          } else {

            htmltools::tags$div(
              style = "text-align: center; margin: 10px 0;",
              htmltools::tags$p(
                style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                paste0(ifelse(metric=="sens","$$ \\hat{C}((\\hat{Se}_{FPR=e})_{1}, (\\hat{Se}_{FPR=e})_{2})","$$ \\hat{C}((\\hat{Sp}_{FNR=e})_{1}, (\\hat{Sp}_{FNR=e})_{2})")," = 0 $$")

              )
            )

          }


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

        htmltools::tagList(

          htmltools::tags$p(
            "Next, compute the variance function of AUC:"
          ),

          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",
            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
              sprintf(
                paste0("$$ \\hat{V}(\\hat{A}) = %.4f \\times e^{-%.3f^2/2} \\times (%.3f \\times%.3f^2+ %.3f +\\dfrac{%.3f \\times %.3f^2+%.3f}{%.3f})= %.3f $$"),
                1/32/pi, a_null, b_null^4+(1+b_null^2)^2, a_null, 2*(1+b_null^2)^2, b_null^4, a_null, 2*b_null^2*(1+b_null^2)^2, R, V0_A1
              )
            )
          ),

          if (paired) {

            htmltools::tagList(

              htmltools::tags$p(
                "Next, compute the covariance function:"
              ),

              htmltools::tags$div(
                style = "text-align: center; margin: 10px 0;",
                htmltools::tags$p(
                  style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                  sprintf(
                    paste0("$$ \\hat{C}_{0}(\\hat{A}_{1}, \\hat{A}_{2})= \\frac{e^{-(\\frac{%.3f}{%.3f}+\\frac{%.3f}{%.3f})}}{2 \\pi \\sqrt{%.3f \\times %.3f}} \\{ (%.3f + \\frac{%.3f}{%.3f}+\\frac{%.3f}{2}) + \\frac{%.3f \\times (%.3f^2+%.3f \\times %.3f^2)}{2\\times %.3f \\times %.3f \\times %.3f} - %.3f \\times [\\frac{%.3f^2}{2 \\times (1+%.3f^2)}+\\frac{%.3f^2}{2 \\times (1+%.3f^2)}] \\} = %.3f $$"),
                    a1_alter^2, (b1_alter^2+1), a2_alter^2, (b2_alter^2+1), (1+b1_alter^2), (1+b2_alter^2), rD, rN*b1_alter*b2_alter, R, rD^2*a1_alter*a2_alter, a1_alter*a2_alter*b1_alter^2*b2_alter^2, rN, R, rD, R, (1+b1_alter^2), (1+b2_alter^2), rD^2*a1_alter*a2_alter, b2_alter, b2_alter, b1_alter, b1_alter, C0_delta_A
                  )
                )
              )
            )

          },

          htmltools::tags$p(
            "Thus, the variance function under the null hypothesis is:"
          ),

          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",
            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
              sprintf(
                paste0("$$ V_{0}(\\hat{A}_{1}-\\hat{A}_{2}) = %.3f + %.3f - 2 \\times %.3f = %.3f $$"),
                V0_A1, V0_A2, C0_delta_A, V0_delta_A
              )
            )
          )
        ),


        htmltools::tags$p(
          "Similarly, compute the variance function term under the alternative hypothesis:"
        ),

        htmltools::tagList(
          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",
            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
              sprintf(
                paste0("$$ V_{A}(\\hat{A}_{1}) = %.3f, V_{A}(\\hat{A}_{2}) = %.3f, C_{A}(\\hat{A}_{1}-\\hat{A}_{2}) = %.3f $$"),
                VA_A1, VA_A2, CA_delta_A
              )
            )
          ),

          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",
            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
              sprintf(
                paste0("$$ V_{A}(\\hat{A}_{1}-\\hat{A}_{2}) = %.3f + %.3f - 2 \\times %.3f = %.3f  $$"),
                VA_A1, VA_A2, CA_delta_A, VA_delta_A
              )
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
                V0_delta_A,
                qnorm(beta),
                VA_delta_A,
                A1_alter,
                A2_alter,
                N
              )
            )

          )

        ),


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
                  paste0("$$ N_{1+} = N_{2+} = n = %d $$"),
                  n_test1_with_condition
                )
              )
            ),

            htmltools::tags$div(
              style = "text-align: center; margin: 10px 0;",
              htmltools::tags$p(
                style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                sprintf(
                  paste0("$$ N_{1-} = N_{2-} = R \\times n = %d $$"),
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
                  paste0("$$ N_{1+} = n = %d, N_{2+} = n = %d$$"),
                  n_test1_with_condition, n_test2_with_condition
                )
              )
            ),

            htmltools::tags$div(
              style = "text-align: center; margin: 10px 0;",
              htmltools::tags$p(
                style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                sprintf(
                  paste0("$$ N_{1-} = R \\times n = %d, N_{2-} = R \\times n = %d$$"),
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

      htmltools::tags$h4(style = "margin-top: 0;", "Participants for test 1 / test 2 with condition"),

      htmltools::tags$p(
        style = "font-size: 24px; font-weight: bold; color: #2E86C1;",
        paste0(n_test1_with_condition," / ", n_test2_with_condition)
      ),

      htmltools::tags$h4(style = "margin-top: 0;", "Participants for test 1 / test 2 without condition"),

      htmltools::tags$p(
        style = "font-size: 24px; font-weight: bold; color: #2E86C1;",
        paste0(n_test1_without_condition," / ", n_test2_without_condition)
      ),



      htmltools::tags$p(
        paste0("This provides ", (1-beta)*100, "% power "),
        paste0(" with ", (1-alpha)*100, "% confidence to detect a difference in area under the ROC curve")
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
        "Sample size calculations for comparing area under the ROC curve"
      ),

      ## 2.Assumptions
      htmltools::tags$p(
        style = "margin: 15px 0px; text-indent: -1em; padding-left: 1.5em;",
        htmltools::tags$strong("2. Assumptions:"),
        htmltools::tags$br(),
        paste0(
          "Assume that ",
          ifelse(dist=="binorm",
                 "the unobserved, underlying test results follow a binormal distribution, but the observed test results, either continuous or ordinal, do not necessarily have a binormal distribution",
                 "the observed test results are on a truly continuous scale and they follow a binormal distribution (or can be transformed to a binormal distribution)"
          )
        )

      ),

      ## 3.Limitations
      htmltools::tags$p(
        style = "margin: 15px 0px; text-indent: -1em; padding-left: 1.5em;",
        htmltools::tags$strong("3. Limitations:"),
        htmltools::tags$br(),
        "Valid on the specific assumptions"
      ),

      ## 4.Reference:
      if(dist == "binorm"){

        htmltools::tags$p(
          style = "margin: 15px 0px; line-height: 1.2;text-indent: -1em; padding-left: 1.5em;",
          htmltools::tags$strong("4. Reference:"),

          htmltools::tags$br(),
          "[1] Obuchowski NA. Computing Sample Size for Receiver Operating Characteristic Studies. ",
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
          "[2] Obuchowski NA, McClish DK. Sample Size Determination for Diagnostic Accuracy Studies Involving Binormal ROC Curve Indices. ",
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
          ),

          htmltools::tags$br(),
          "[3] Rockette HE, Campbell WL, Britton CA, Holbert JM, King JL, Gur D. Empiric Assessment of Parameters That Affect the Design of Multireader Receiver Operating Characteristic Studies. ",
          htmltools::tags$em("Academic Radiology."),
          " 1999;6(12):723-729. ",
          htmltools::tags$a(
            href = "doi:https://doi.org/10.1016/s1076-6332(99)80468-1",
            target = "_blank",
            style = "color: #1a5fb4; text-decoration: none; border-bottom: 1px solid #c1dbf7;
             transition: border-color 0.3s, color 0.3s;",
            onmouseover = "this.style.borderBottomColor='#1a5fb4'; this.style.color='#0d4e9b';",
            onmouseout = "this.style.borderBottomColor='#c1dbf7'; this.style.color='#1a5fb4';",
            "doi:https://doi.org/10.1016/s1076-6332(99)80468-1"
          )

        )
      } else if (dist == "binorm"){




      }
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
          "to compare the primary diagnostic performance metrics (e.g., AUC and partial AUC) in ",
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
        parameters = list(A_null = A_null, b_null = b_null, A1_alter = A1_alter, b1_alter = b1_alter, A2_alter = A2_alter, b2_alter = b2_alter, e1 = e1, e2 = e2, alpha = alpha, beta = beta, rN = rN, rD = rD, R = R, paired = paired),
        html_report = html_report
      ),
      class = "diag_sample_size_html"
    )
  )


}
