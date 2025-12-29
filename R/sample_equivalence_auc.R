#' Sample Size Calculator for Diagnostic Tests
#'
#' Calculate the required sample size for testing whether the Area Under the ROC Curve (AUC) of a new test is equivalent to an existing test.
#'
#' @param As (numeric): The expected sensitivity (or specificity) of an existing test (standard test).
#' @param Ae (numeric): The expected sensitivity (or specificity) of a new test (experimental test).
#' @param prob_t1pos_give_t2pos (numeric): Conditional probability of the existing test being positive when the new test being positive.
#' @param alpha (numeric): Significance level (default: 0.05).
#' @param beta (numeric): Type II error rate.
#' @param delta (numeric):  The prespecified upper clinical limits for equivalence.
#' @param paired (logical): Whether the design is paired (`TRUE`) or unpaired (`FALSE`).
#' @param dist (character): The assumption for choosing the variance function (default: "any"):
#'   - `"any"`: Applicable for all distributions.
#'   - `"binorm"`: Assume that the test results have an underlying bivariate binormal distribution.
#' @return An object of class "diag_sample_size_html" containing:
#'   - `sample_size` (list): Required sample size.
#'   - `parameters` (list): Input parameters.
#'   - `html_report` (html): Sample size calculation report.
#' @examples
#' n <- sample_equivalence_auc(As=0.9, bs=1, Ae=0.88, be=1, rN=0.5, rD=0.5, delta=0.05, alpha=0.05, beta=0.2, R=1, paired=TRUE, dist="binorm")
#' n <- sample_equivalence_auc(As=0.9, bs=1, Ae=0.90, be=1, rN=0.5, rD=0.5, delta=0.05, alpha=0.05, beta=0.2, R=1, paired=TRUE, dist="binorm")
#' n <- sample_equivalence_auc(As=0.9, bs=1, Ae=0.92, be=1, rN=0.5, rD=0.5, delta=0.05, alpha=0.05, beta=0.2, R=1, paired=TRUE, dist="binorm")
#' print(n)
#' @export
sample_equivalence_auc <- function(As, bs=1, Ae, be=1, r=NULL, rN=NULL, rD=NULL, delta, alpha=0.05, beta, R=1, paired=TRUE, dist="any") {

  #  (1) ----- validate inputs

  if (!is.numeric(As) || length(As) != 1 || As <= 0 ||As >= 1) {
    stop("'As' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(Ae) || length(Ae) != 1 || Ae <= 0 || Ae >= 1) {
    stop("'Ae' must be a single numeric value in (0, 1).")
  }

  if(dist == "any") {

    if (!is.numeric(r) || length(r) != 1 || r < -1 || r > 1) {
      stop("'r' must be a single numeric value in [-1, 1].")
    }

  }
  else if (dist == "binorm") {

    if (!is.numeric(rN) || length(rN) != 1 || rN < -1 || rN > 1) {
      stop("'rN' must be a single numeric value in [-1, 1].")
    }

    if (!is.numeric(rD) || length(rD) != 1 || rD < -1 || rD > 1) {
      stop("'rD' must be a single numeric value in [-1, 1].")
    }

    if (!is.numeric(bs) || length(bs) != 1 || bs <= 0) {
      stop("'bs' must be a single positive numeric value.")
    }

    if (!is.numeric(be) || length(be) != 1 || be <= 0) {
      stop("'be' must be a single positive numeric value.")
    }

  } else {
    stop("'dist' must be either 'any' or 'binorm'.")

  }

  if (!is.numeric(alpha) || length(alpha) != 1 || alpha <= 0 || alpha >= 1) {
    stop("'alpha' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(beta) || length(beta) != 1 || beta <= 0 || beta >= 1) {
    stop("'beta' must be a single numeric value in (0, 1).")
  }

  if (!is.numeric(delta) || length(delta) != 1 || delta <= 0 || delta >= 1) {
    stop("'delta' must be a single numeric value in (0, 1).")
  }

  if (!is.logical(paired) || length(paired) != 1) {
    stop("'paired' must be a single logical value (TRUE or FALSE).")
  }


  #  (2) -----  calculation modules
  message(paste0("Null hypothesis: θs - θe ≤ ", -delta, " or θs - θe ≥ ", delta))
  message(paste0("Alternative hypothesis: ", -delta, " < θs - θe <", delta))

  ## define correlation coefficients for paired and unpaired designs

  if(paired){

    message("For paired design, r = rD = rN = 0.5 is a typical value for the correlation between tests")

    if(dist == "any"){

      message("If not specified, r defaults to 0.5")

      rN <- NA
      rD <- NA

    } else if (dist == "binorm") {

      message("If not specified, rD and rN default to 0.5")

      r <- NA

    }

  } else {

    message("For unpaired design, r = rD = rN = 0")
    r <- 0
    rN <- rD <- r

  }


  ## calculate variance function
  as <- get_binorm_params(As, bs)$a
  ae <- get_binorm_params(Ae, be)$a

  if(dist == "any") {

    VA_delta_A <- As*(1-As)+Ae*(1-Ae)-2*r*sqrt(As*(1-As)*Ae*(1-Ae))
    var_function <- VA_delta_A

  } else if (dist == "binorm") {

    VA <- function(a, b) {

      # V <- (0.0099) * exp(-a^2 / 2) * ((5 * a^2 + 8) + (a^2 + 8) / R) # if b = 1
      expr1 <- exp(- a^2 / 2 / (1 + b^2))
      expr2 <- 1 + b^2

      f <- expr1 * (2 * pi * expr2)^(-1/2)
      g <- -expr1 * (a*b) * (2 * pi * expr2^3)^(-1/2)

      V <- f^2 * (1 + b^2 / R + a^2 / 2) + g^2 * (b^2 * (1 + R)/ 2 / R)
      return(V)

    }

    VA_As <- VA(as, bs)
    VA_Ae <- VA(ae, be)

    C_delta_A <- function(a1, a2, b1, b2, rD, rN, R) {

      expr11 <- exp(- a1^2 / 2 / (1 + b1^2))
      expr12 <- 1 + b1^2
      f1 <- expr11 * (2 * pi * expr12)^(-1/2)
      g1 <- -expr11 * (a1*b1) * (2 * pi * expr12^3)^(-1/2)

      expr21 <- exp(- a2^2 / 2 / (1 + b2^2))
      expr22 <- 1 + b2^2
      f2 <- expr21 * (2 * pi * expr22)^(-1/2)
      g2 <- -expr21 * (a2*b2) * (2 * pi * expr22^3)^(-1/2)

      C <- f1*f2*(rD + rN*b1*b2/R + rD^2*a1*a2/2) + g1*g2*(b1*b2*(rN^2+R*rD^2)/2/R) + f1*g2*rD^2*a1*b2/2 + f2*g1*rD^2*a2*b1/2
      # term1 <- exp(-(a1^2 + a2^2) / 4) / 12.5664 * (rD + rN / R + rD^2* a1 * a2 / 2)
      # term2 <- exp(-(a1^2 + a2^2) / 4) / 50.2655 * a1 * a2 * (rN^2 + R* rD^2) / 2 / R
      # term3 <- - exp(-(a1^2+a2^2) / 4) / 25.1327 * rD^2 * a1 * a2
      # C <- term1 + term2 + term3
      return(C)
    }

    CA_delta_A <- ifelse(paired, C_delta_A(as, ae, bs, be, rD, rN, R), 0)

    VA_delta_A <- VA_As + VA_Ae - 2 * CA_delta_A
    var_function <- VA_delta_A

  }


  ## calculate delta
  if(As - Ae > 0) {

    delta1 <- abs(delta - (As - Ae))

  } else if (As - Ae < 0) {

    delta1 <- abs(delta + (As - Ae))

  } else if (As - Ae == 0){

    beta = beta / 2
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
        htmltools::tags$strong(style = paste0("color:", ifelse(as > 0.5, "#E74C3C", "#27AE60"), ";"),"area under the ROC curve"),
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
                      "Binominal Parameter the Existing Test",
                      "Binominal Parameter the New Test",
                      "Correlation Between Test",
                      "Correlation for Patients with Condition",
                      "Correlation for Patients without Condition",
                      "Confidence Level",
                      "Statistical Power",
                      "Equivalence Margin",
                      "Group Allocation"

        ),
        Symbol    = c("As",
                      "Ae",
                      "(as, bs)",
                      "(ae, be)",
                      "r",
                      "rN",
                      "rD",
                      "1 - α",
                      "1 - β",
                      "(ΔL, ΔU)",
                      "R"
        ),
        Value     = c(As,
                      Ae,
                      paste0("(", round(as,3), ", ", round(bs,3), ")"),
                      paste0("(", round(ae,3), ", ", round(be,3), ")"),
                      ifelse(is.na(r), "Not specified", r),
                      ifelse(is.na(rN), "Not specified", rN),
                      ifelse(is.na(rD), "Not specified", rD),
                      paste0((1-alpha)*100, "%"),
                      paste0((1-beta)*100, "%"),
                      paste0("(",-delta,", ",delta,")"),
                      R
        ),
        Notes     = c("The expected AUC of the existing test",
                      "The expected AUC of the new test",
                      "a = (μs_with - μs_without)/σs_with; b = σs_without/σs_with, to determine the shape of ROC of the existing test",
                      "a = (μe_with - μe_without)/σe_with; b = σe_without/σe_with, to determine the shape of ROC of the new test",
                      "The correlation between the tests",
                      "The correlation of the underlying bivariate binormal distribution for patients with the condition",
                      "The correlation of the underlying bivariate binormal distribution for patients without the condition",
                      "1 - Type I error rate",
                      "1 - Type II error rate",
                      "The prespecified lower and upper clinical limits for equivalence",
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
          "For the area under the ROC curve (AUC) equivalence testing, the null and alternative hypotheses are formulated as:"
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

            "$$ H_{0}: A_{S} - A_{E} \\leq \\Delta_{L} \\quad\\text{or}\\quad A_{S} - A_{E} \\geq \\Delta_{U} $$",
            "$$ H_{1}: \\Delta_{L} < A_{S} - A_{E} < \\Delta_{U} $$"
          )
        ),

        htmltools::tags$p(
          "The required sample size is calculated using the formula:"
        ),

        if (As - Ae == 0) {

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

              "$$ n = \\dfrac{(Z_{\\alpha}+Z_{\\beta/2})^2 V_{A}(\\hat{A}_{S}-\\hat{A}_{E})}{(Δ_{U})^2} $$"
            )
          )

        } else if(As - Ae != 0) {

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

              "$$ n = \\dfrac{(Z_{\\alpha}+Z_{\\beta})^2 V_{A}(\\hat{A}_{S}-\\hat{A}_{E})}{(Δ_{U}-|A_{S}-A_{E}|)^2} $$"
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

        if(dist == "any") {

          htmltools::tagList(

            htmltools::tags$p(
              "The variance function of the estimated difference in AUC under the alternative hypothesis is calculated through:"
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

                "$$ V_{A}(\\hat{A}_{S} - \\hat{A}_{E}) = A_{S}(1-A_{S})+A_{E}(1-A_{E}) - 2r\\sqrt{A_{S}(1-A_{S})A_{E}(1-A_{E})} $$"
              )
            )
          )

        } else if(dist == "binorm") {

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

                "$$ V(\\hat{A}_{S}-\\hat{A}_{E}) = V(\\hat{A}_{S}) + V(\\hat{A}_{E}) - 2C(\\hat{A}_{S},\\hat{A}_{E}) $$"
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
                    paste0("$$ \\hat{C}(\\hat{A}_{S}, \\hat{A}_{E})= \\frac{e^{-[\\frac{a_{S}^2}{2(b_{S}^2+1)}+\\frac{a_{E}^2}{2(b_{E}^2+1)}]}}{2\\pi\\sqrt{(1+b_{S}^2)(1+b_{E}^2)}} \\{ (r_{D} + \\frac{r_{N}b_{S}b_{E}}{R}+\\frac{r_{D}^2a_{S}a_{E}}{2}) + \\frac{a_{S}a_{E}b_{S}^2b_{E}^2(r_{N}^2+Rr_{D}^2)}{2R(1+b_{S}^2)(1+b_{E}^2)} - r_{D}^2a_{S}a_{E}[\\frac{b_{E}^2}{2(1+b_{E}^2)}+\\frac{b_{S}^2}{2(1+b_{S}^2)}] \\} $$")

                  )
                ),

                htmltools::tags$p("Specifically, when b", htmltools::tags$sub("S"), " = b", htmltools::tags$sub("E"), " = 1"),

                htmltools::tags$div(
                  style = "text-align: center; margin: 10px 0;",
                  htmltools::tags$p(
                    style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                    paste0("$$ \\hat{C}(\\hat{A}_{S}, \\hat{A}_{E})= e^{-(a_{S}^2+a_{E}^2)/4} [\\frac{r_{D} + \\frac{r_{N}}{R}+\\frac{r_{D}^2a_{S}a_{E}}{2}}{12.5664} + \\frac{a_{S}a_{E}(\\frac{r_{N}^2}{R}+r_{D}^2)}{100.531} - \\frac{r_{D}^2a_{S}a_{E}}{25.1327}] $$")

                  )
                )
              )

            } else {

              htmltools::tags$div(
                style = "text-align: center; margin: 10px 0;",
                htmltools::tags$p(
                  style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                  paste0("$$ \\hat{C}(\\hat{A}_{S}, \\hat{E}_{2})= 0 $$")

                )
              )

            }


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


        if(dist == "any") {

          htmltools::tagList(
            htmltools::tags$p(
              "Next, compute the variance function term under the alternative hypothesis:"
            ),

            htmltools::tags$div(
              style = "text-align: center; margin: 10px 0;",
              htmltools::tags$p(
                style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                sprintf(
                  paste0("$$ V_{A}(\\hat{A}_{S} - \\hat{A}_{E}) = %.3f \\times (1-%.3f) + %.3f \\times (1-%.3f) - 2 \\times %.3f\\sqrt{%.3f \\times (1-%.3f) \\times %.3f \\times (1-%.3f)} = %.3f $$"),
                  As, As, Ae, Ae, r, As, As, Ae, Ae, VA_delta_A
                )
              )
            )

          )

        } else if(dist == "binorm"){

          htmltools::tagList(

            htmltools::tags$p(
              "Next, compute the variance function of ",
              "A", htmltools::tags$sub("s"), " and ",
              "A", htmltools::tags$sub("e"), ":"
            ),

            htmltools::tags$div(
              style = "text-align: center; margin: 10px 0;",
              htmltools::tags$p(
                style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                sprintf(
                  paste0("$$ \\hat{V}(\\hat{A}_{S}) = %.4f \\times e^{-%.3f^2/2} \\times (%.3f \\times%.3f^2+ %.3f +\\dfrac{%.3f \\times %.3f^2+%.3f}{%.3f})= %.3f $$"),
                  1/32/pi, as, bs^4+(1+bs^2)^2, as, 2*(1+bs^2)^2, bs^4, as, 2*bs^2*(1+bs^2)^2, R, VA_As
                )
              )
            ),

            htmltools::tags$div(
              style = "text-align: center; margin: 10px 0;",
              htmltools::tags$p(
                style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                sprintf(
                  paste0("$$ \\hat{V}(\\hat{A}_{E}) = %.4f \\times e^{-%.3f^2/2} \\times (%.3f \\times%.3f^2+ %.3f +\\dfrac{%.3f \\times %.3f^2+%.3f}{%.3f})= %.3f $$"),
                  1/32/pi, ae, be^4+(1+be^2)^2, ae, 2*(1+be^2)^2, be^4, ae, 2*be^2*(1+be^2)^2, R, VA_Ae
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
                      paste0("$$ \\hat{C}_{A}(\\hat{A}_{S}, \\hat{A}_{E})= \\frac{e^{-(\\frac{%.3f}{%.3f}+\\frac{%.3f}{%.3f})}}{2 \\pi \\sqrt{%.3f \\times %.3f}} \\{ (%.3f + \\frac{%.3f}{%.3f}+\\frac{%.3f}{2}) + \\frac{%.3f \\times (%.3f^2+%.3f \\times %.3f^2)}{2\\times %.3f \\times %.3f \\times %.3f} - %.3f \\times [\\frac{%.3f^2}{2 \\times (1+%.3f^2)}+\\frac{%.3f^2}{2 \\times (1+%.3f^2)}] \\} = %.3f $$"),
                      as^2, (bs^2+1), ae^2, (be^2+1), (1+bs^2), (1+be^2), rD, rN*bs*be, R, rD^2*as*ae, as*ae*bs^2*be^2, rN, R, rD, R, (1+bs^2), (1+be^2), rD^2*as*ae, be, be, bs, bs, CA_delta_A
                    )
                  )
                )
              )

            },

            htmltools::tags$p(
              "Thus, the variance function under the alternative hypothesis is:"
            ),

            htmltools::tags$div(
              style = "text-align: center; margin: 10px 0;",
              htmltools::tags$p(
                style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",
                sprintf(
                  paste0("$$ V_{A}(\\hat{A}_{S}-\\hat{A}_{E}) = %.3f + %.3f - 2 \\times %.3f = %.3f $$"),
                  VA_As, VA_Ae, CA_delta_A, VA_delta_A
                )
              )
            )
          )

        },


        htmltools::tags$p(
          "Now apply the full formula:"
        ),

        if(As - Ae == 0) {

          htmltools::tags$div(
            style = "text-align: center; margin: 15px 0;",

            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",

              htmltools::HTML(
                sprintf(
                  "$$n = \\dfrac{[(%.3f) + (%.3f)]^2 \\times %.3f}{(%.3f)^2} ≈ %d$$",
                  qnorm(alpha),
                  qnorm(beta),
                  VA_delta_A,
                  delta,
                  N
                )
              )

            )

          )
        } else if (As - Ae != 0) {

          htmltools::tags$div(
            style = "text-align: center; margin: 15px 0;",

            htmltools::tags$p(
              style = "font-family: 'Cambria Math', serif; font-size: 1.1em;",

              htmltools::HTML(
                sprintf(
                  "$$n = \\dfrac{[(%.3f) + (%.3f)]^2 \\times %.3f}{(%.3f - |%.3f-%.3f|)^2} ≈ %d$$",
                  qnorm(alpha),
                  qnorm(beta),
                  VA_delta_A,
                  delta,
                  As,
                  Ae,
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
        paste0(" with ", (1-alpha)*100, "% confidence to test equivalence in area under the ROC curve")
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
          ),

          htmltools::tags$br(),
          "[4] Dongsheng T. Two one-sided Tests Procedures in Establishing Therapeutic Equivalence with Binary Clinical endpoints: Fixed Sample Performances and Sample Size Derformination. ",
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
          "to demonstrate that the difference in primary diagnostic performance metrics (e.g., AUC and partial AUC) falls within the equivalence margin in ",
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
        parameters = list(As = As, bs = bs, Ae = Ae, be=be, alpha = alpha, beta = beta, r = r, rN = rN, rD = rD, R = R, paired = paired),
        html_report = html_report
      ),
      class = "diag_sample_size_html"
    )
  )


}
