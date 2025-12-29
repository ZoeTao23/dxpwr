#' Sample Size Calculator for Diagnostic Tests
#'
#' Calculate the required sample size for comparing partial Area Under the ROC Curve (pAUC) between two diagnostic tests.
#'
#' @param A_null (numeric): AUC under the null hypothesis.
#' @param b_null (numeric): Ratio of the standard deviations of the distributions of test results for patients without versus with the condition, σ_without/σ_with (default: 1).
#' @param A1_alter (numeric): AUC of test 1 under the alternative hypothesis.
#' @param b1_alter (numeric): Ratio of the standard deviations of the distributions of results for test 1 for patients without versus with the condition, σ1_without/σ1_with (default: 1).
#' @param A2_alter (numeric): AUC of test 2 under the alternative hypothesis.
#' @param b2_alter (numeric): Ratio of the standard deviations of the distributions of results for test 2 for patients without versus with the condition, σ2_without/σ2_with (default: 1).
#' @param e1 (numeric): Lower bound of the false positive rate (FPR) range.
#' @param e2 (numeric): Upper bound of the false positive rate (FPR) range.
#' @param rD (numeric): The correlation of the underlying bivariate binormal distribution for patients with the condition (default: 0.5).
#' @param rN (numeric): The correlation of the underlying bivariate binormal distribution for patients without the condition (default: 0.5).
#' @param alpha (numeric): Significance level (default: 0.05).
#' @param beta (numeric): Type II error rate.
#' @param paired (logical). Whether the design is paired (`TRUE`) or unpaired (`FALSE`).
#' @param R (numeric): Ratio of patients with and without the condition (default: 1).
#' @return An object of class "diag_sample_size_html" containing:
#'   - `sample_size` (list): Required sample size.
#'   - `parameters` (list): Input parameters.
#'   - `html_report` (html): Sample size calculation report.
#' @examples
#' n <- sample_compare_partial_roc(A_null= 0.8, b_null=1, A1_alter=0.8, b1_alter=1, A2_alter=0.9, b2_alter=1, e1=0, e2=0.2, rN=0, rD=0, alpha = 0.05, beta=0.2, paired=FALSE)
#' print(n)
#' @export
sample_compare_partial_roc <- function(A_null, b_null=1, A1_alter, b1_alter=1, A2_alter, b2_alter=1, e1, e2, rN=0.5, rD=0.5, alpha=0.05, beta, R=1, paired=FALSE) {


  #  (1) ----- validate inputs
  if (!is.numeric(A_null) || length(A_null) != 1 || A_null <= 0 ||A_null >= 1) {
    stop("'A_null' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(A1_alter) || length(A1_alter) != 1 || A1_alter <= 0 ||A1_alter >= 1) {
    stop("'A1_alter' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(A2_alter) || length(A2_alter) != 1 || A2_alter <= 0 ||A2_alter >= 1) {
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

  if (!is.numeric(e1) || length(e1) != 1 || e1 < 0 || e1 >= 1) {
    stop("'e1' must be a single numeric value in [0, 1).")
  }

  if (!is.numeric(e2) || length(e2) != 1 || e2 <= 0 || e2 > 1) {
    stop("'e2' must be a single numeric value in (0, 1].")
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

  if (!is.numeric(beta) || length(beta) != 1 || beta <= 0 || beta >= 1) {
    stop("'beta' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(R) || length(R) != 1 || R <= 0) {
    stop("'R' must be a single positive numeric value.")
  }

  if (!is.logical(paired) || length(paired) != 1) {
    stop("'paired' must be a single logical value (TRUE or FALSE).")
  }

  #  (2) ----- calculation modules

  message("Null hypothesis: AUC1 = AUC2 when fdr is between (",e1,",",e2,")")
  message("Alternative hypothesis: AUC1 = AUC2 = ", A_null, " when fdr is between (",e1,",",e2,")")

  ## define the correlation
  if(paired){

    message("For paired design, rD = rN = 0.5 is a typical value for the correlation between tests")
    message("If not specified, rD and rN default to 0.5")

  } else if(!paired) {

    message("For unpaired design, rD = rN = 0")
    r <- 0
    rN <- rD <- r
  }

  ## calculate variance function
  f <- function(a, b, e1, e2){

    e_prime <- (qnorm(c(e1,e2)) + a * b / (1 + b^2 )) * sqrt(1 + b^2)

    expr1 <- exp(-a^2 / 2 / (1 + b^2))
    expr2 <- 1 + b^2
    expr3 <- pnorm(e_prime[2]) - pnorm(e_prime[1])

    f_value <- expr1 * (2 * pi * expr2) ^ (-1/2) * expr3
    return(f_value)

  }


  g <- function(a, b, e1, e2){

    e_prime <- (qnorm(c(e1,e2)) + a * b / (1 + b^2)) * sqrt(1 + b^2)
    e_double_prime <- (e_prime)^2 / 2

    expr1 <- exp(-a^2 / 2 / (1 + b^2))
    expr2 <- 1 + b^2
    expr3 <- pnorm(e_prime[2]) - pnorm(e_prime[1])
    expr4 <- exp(-e_double_prime[1]) - exp(-e_double_prime[2])

    g_value <- expr1 * (2 * pi * expr2)^(-1) * expr4 - a * b* expr1 * (2 * pi * expr2^3)^(-1/2) * expr3
    return(g_value)

  }


  VA <- function(a, b, e1, e2) {

    f <- f(a, b, e1, e2)
    g <- g(a, b, e1, e2)
    V <- f^2 * (1 + b^2 / R + a^2 / 2) + g^2 * (b^2 * (1+R) / 2 / R)

    return(V)
  }

  C_delta_A <- function(a1, a2, b1, b2, e1, e2, rN, rD) {

    f1 <- f(a1, b1, e1, e2)
    f2 <- f(a2, b2, e1, e2)
    g1 <- g(a1, b1, e1, e2)
    g2 <- g(a2, b2, e1, e2)

    term1 <- f1 * f2 * (rD + rN * b1 * b2 / R + rD^2 * a1 * a2 / 2)
    term2 <- g1 * g2 * b1 * b2 * (rN^2 + R * rD^2) / 2 / R
    term3 <- f1 * g2 * rD^2 * a1 * b2 / 2
    term4 <- f2 * g1 * rD^2 * a2 * b1

    C <- term1 + term2 + term3 + term4
    return(C)

  }

  a_null <- get_binorm_params(A_null, b_null)$a
  a1_alter <- get_binorm_params(A1_alter, b1_alter)$a
  a2_alter <- get_binorm_params(A2_alter, b2_alter)$a

  f_null <- f(a_null, b_null, e1, e2)
  g_null <- g(a_null, b_null, e1, e2)
  f1 <- f(a1_alter, b1_alter, e1, e2)
  g1 <- g(a1_alter, b1_alter, e1, e2)
  f2 <- f(a2_alter, b2_alter, e1, e2)
  g2 <- g(a2_alter, b2_alter, e1, e2)


  ## calculate variance function
  V0_A1 <- V0_A2 <- VA(a_null, b_null, e1, e2)
  C0_delta_A <- C_delta_A(a_null, a_null, b_null, b_null, e1, e2, rN, rD)
  V0_delta_A <- V0_A1 + V0_A2 - 2 * C0_delta_A

  VA_A1 <- VA(a1_alter, b1_alter, e1, e2)
  VA_A2 <- VA(a2_alter, b2_alter, e1, e2)
  CA_delta_A <- C_delta_A(a1_alter, a2_alter, b1_alter, b2_alter, e1, e2, rN, rD)
  VA_delta_A <- VA_A1 + VA_A2 - 2 * CA_delta_A

  var_function <- c(V0_delta_A, VA_delta_A)

  ## calculate delta
  fpr <- function(a, b, x) {
    e <- pnorm(a + b * qnorm(x))
    return(e)
  }

  delta1 <- abs(integrate(function(x) fpr(a1_alter, b1_alter, x), e1, e2)$value -
                  integrate(function(x) fpr(a2_alter, b2_alter, x), e1, e2)$value)


  ## calculate sample size
  N <- get_diag_sample(
    var_function = var_function,
    alpha = alpha,
    beta = beta,
    delta = delta1,
    test_type = "two_diagnostic"
  )


  ## calculate sample size with and without condition
  if(paired){
    n_test1_with_condition <- ceiling(N)
    n_test1_without_condition <- ceiling(N*R)
    n_test1 <- n_test1_with_condition + n_test1_without_condition

    n_test2_with_condition <- ceiling(N)
    n_test2_without_condition <- ceiling(N*R)
    n_test2 <- n_test2_with_condition + n_test2_without_condition

    n_total <- n_test1

  } else {
    n_test1_with_condition <- ceiling(N)
    n_test1_without_condition <- ceiling(N*R)
    n_test1 <- n_test1_with_condition + n_test1_without_condition

    n_test2_with_condition <- ceiling(N)
    n_test2_without_condition <- ceiling(N*R)
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
        "To compare the ",
        htmltools::tags$strong(style = paste0("color:", ifelse(a1_alter > 0.5, "#E74C3C", "#27AE60"), ";"),"patial area under the ROC curve"),
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
                      "FPR Range",
                      "Correlation for Patients with Condition",
                      "Correlation for Patients without Condition",
                      "Confidence Level",
                      "Statistical Power",
                      "Difference in AUC",
                      "Group Allocation"

        ),
        Symbol    = c("A_null",
                      "A1_alter",
                      "A2_alter",
                      "(a_null, b_null)",
                      "(a1_alter, b1_alter)",
                      "(a2_alter, b2_alter)",
                      "(e1, e2)",
                      "rN",
                      "rD",
                      "1 - α",
                      "1 - β",
                      "Δ1",
                      "R"
        ),
        Value     = c(ifelse(is.null(A_null),"Not specified", A_null),
                      A1_alter,
                      A2_alter,
                      ifelse(is.null(A_null),"Not specified", paste0("(", round(a_null,3), ", ", round(b_null,3), ")")),
                      paste0("(", round(a1_alter,3), ", ", round(b1_alter,3), ")"),
                      paste0("(", round(a2_alter,3), ", ", round(b2_alter,3), ")"),
                      paste0("(", e1, ", ", e2, ")"),
                      ifelse(is.na(rN), "Not specified", rN),
                      ifelse(is.na(rD), "Not specified", rD),
                      paste0((1-alpha)*100, "%"),
                      paste0((1-beta)*100, "%"),
                      abs(A1_alter-A2_alter),
                      R
        ),
        Notes     = c("The expected AUC in null hypothesis",
                      "The expected AUC of test 1 in alternative hypothesis",
                      "The expected AUC of test 2 in alternative hypothesis",
                      "a = (μ0_with - μ0_without)/σ0_with; b = σ0_without/σ0_with, to determine the shape of ROC in null hypothesis",
                      "a = (μ1_with - μ1_without)/σ1_with; b = σ1_without/σ1_with, to determine the shape of ROC of test 1 in alternative hypothesis",
                      "a = (μ2_with - μ2_without)/σ2_with; b = σ2_without/σ2_with, to determine the shape of ROC of test 2 in alternative hypothesis",
                      "e1 = Lower bound of FPR, e2 = Upper bound of FPR",
                      "The correlation of the underlying bivariate binormal distribution for patients with the condition",
                      "The correlation of the underlying bivariate binormal distribution for patients without the condition",
                      "1 - Type I error rate",
                      "1 - Type II error rate",
                      "Difference between AUC under the alternative hypothesis",
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

        htmltools::tags$p(
          "For comparing patrial area under the ROC curve (AUC) with a given precision, the null and alternative hypotheses are:"
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

            "$$ H_{0}: AUC_{1,e_{1} \\leq FPR \\leq e_{2}} = AUC_{2,e_{1} \\leq FPR \\leq e_{2}}() $$",
            "$$ H_{1}: AUC_{1,e_{1} \\leq FPR \\leq e_{2}}() ≠ AUC_{2,e_{1} \\leq FPR \\leq e_{2}} $$"
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

          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",
            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
              paste0("$$ \\hat{V}(\\hat{A})= \\frac{e^{-a^2/2}}{4\\pi(1+b^2)^3} \\times\\{a^2[b^4+(1+b^2)^2]+2(1+b^2)^2+\\frac{a^2b^4+2b^2(1+b^2)^2}{R}\\} $$")

            )
          ),

          htmltools::tags$p(
            "Specifically, when b = 1"
          ),

          htmltools::tags$div(
            style = "text-align: center; margin: 10px 0;",
            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
              paste0("$$ \\hat{V}(\\hat{A})=0.0099 \\times e^{-a^2/2} \\times(5a^2+8+\\dfrac{a^2+8}{R}) $$")

            )
          ),

          htmltools::tags$p(
            paste0(ifelse(paired, "And for paired design, ", "And for unpaired design, "), "the covariance function is:")
          ),

          if(paired) {

            htmltools::tagList(
              htmltools::tags$div(
                style = "text-align: center; margin: 10px 0;",
                htmltools::tags$p(
                  style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                  paste0("$$ \\hat{C}(\\hat{A}_{1}, \\hat{A}_{2})= \\frac{e^{-[\\frac{a_{1}^2}{2(b_{1}^2+1)}+\\frac{a_{2}^2}{2(b_{2}^2+1)}]}}{2\\pi\\sqrt{(1+b_{1}^2)(1+b_{2}^2)}} \\{ (r_{D} + \\frac{r_{N}b_{1}b_{2}}{R}+\\frac{r_{D}^2a_{1}a_{2}}{2}) + \\frac{a_{1}a_{2}b_{1}^2b_{2}^2(r_{N}^2+Rr_{D}^2)}{2R(1+b_{1}^2)(1+b_{2}^2)} - r_{D}^2a_{1}a_{2}[\\frac{b_{2}^2}{2(1+b_{2}^2)}+\\frac{b_{1}^2}{2(1+b_{1}^2)}] \\} $$")

                )
              ),

              htmltools::tags$p(
                "Specifically, when b1  = b2 = 1"
              ),

              htmltools::tags$div(
                style = "text-align: center; margin: 10px 0;",
                htmltools::tags$p(
                  style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                  paste0("$$ \\hat{C}(\\hat{A}_{1}, \\hat{A}_{2})= e^{-(a_{1}^2+a_{2}^2)/4} [\\frac{r_{D} + \\frac{r_{N}}{R}+\\frac{r_{D}^2a_{1}a_{2}}{2}}{12.5664} + \\frac{a_{1}a_{2}(\\frac{r_{N}^2}{R}+r_{D}^2)}{100.531} - \\frac{r_{D}^2a_{1}a_{2}}{25.1327}] $$")

                )
              )
            )

          } else {

            htmltools::tags$div(
              style = "text-align: center; margin: 10px 0;",
              htmltools::tags$p(
                style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                paste0("$$ \\hat{C}(\\hat{A}_{1}, \\hat{A}_{2})= 0 $$")

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
              "$$z_{\\alpha/2} = z_{%.2f} = %.3f $$",
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
                  paste0("$$ N_{1-} = N_{2-} = n = %d $$"),
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
                  paste0("$$ N_{1-} = Rn = %d, N_{2-} = Rn = %d$$"),
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

      if(dist == "any") {
        htmltools::tags$p(
          style = "margin: 15px 0px; text-indent: -1em; padding-left: 1.5em;",
          htmltools::tags$strong("2. Assumptions:"),
          htmltools::tags$br(),
          "No specific restrictions on the distributions of the test results"

        )
      } else if(dist == "binorm") {

        htmltools::tags$p(
          style = "margin: 15px 0px; text-indent: -1em; padding-left: 1.5em;",
          htmltools::tags$strong("2. Assumptions:"),
          htmltools::tags$br(),
          "Assume that the test results have an underlying bivariate binormal distribution"

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
      if (dist == "any") {

        htmltools::tags$p(
          style = "margin: 15px 0px; line-height: 1.2;text-indent: -1em; padding-left: 1.5em;",
          htmltools::tags$strong("4. Reference:"),

          htmltools::tags$br(),
          "[1] Blume JD. Bounding Sample Size Projections for the Area under a ROC Curve. ",
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
          ),

          htmltools::tags$br(),
          "[2] Rockette HE, Campbell WL, Britton CA, Holbert JM, King JL, Gur D. Empiric Assessment of Parameters That Affect the Design of Multireader Receiver Operating Characteristic Studies. ",
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

      } else if(dist == "binorm"){

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
