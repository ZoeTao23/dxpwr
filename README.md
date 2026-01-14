# dxpwr：An Integrated R Toolkit for Sample Size Determination in Diagnostic Accuracy Studies
![Supported Designs Overview](workflow.png)
## 📋 Overview
`dxpwr` is an R package designed to provide a structured, transparent framework for sample size determination in diagnostic accuracy studies. It supports a wide range of study designs—from single‑test evaluation to comparative studies (difference‑based, equivalence, non‑inferiority)—with built‑in adjustments for practical factors like prospective design and clustered data.
A companion web‑based calculator offers the main functionality through an intuitive, point‑and‑click interface, making rigorous sample‑size planning accessible to users without programming experience.

## ✨ Key Features
1. **Guided, transparent workflow** – Maps your study question (single‑test or comparative, paired or unpaired) directly to the appropriate statistical estimator.
2. **Comprehensive design support** – Covers sensitivity, specificity, AUC, partial AUC, sensitivity at a fixed false‑positive rate, specificity at a fixed false‑negative rate, and more.
3. **Practical adjustments** – Automatically converts case/control numbers into total recruitment size for prospective studies, and inflates sample size for clustered data.
4. **Fully reproducible reporting** – For every calculation, `dxpwr` generates a detailed report that includes the statistical formula used, a breakdown of each input parameter, and a step‑by‑step computational trace.
5. **Dual interface** – Use the R package for programmatic, reproducible workflows, or the web app for interactive, guided planning.

## 🧪 Supported Designs 
![Supported Designs Overview](user_guide.png)
Both **paired** (same participants receive both tests) and **unpaired** (different participants for each test) designs are supported.

## 📦 Installation
```r
install.packages("devtools")
devtools::install_github("ZoeTao23/dxpwr")
```

## 🚀 Quick Start
```r
library(dxpwr)

# Example 1: Calculate sample size for a single-test study

# Number of sample needed when estimating sensitivity of a new test
report <- sample_estimate_sesp(
  theta=0.8,    # Pre-specified sensitivity (or specificity)
  alpha=0.05,   # Type I error rate
  beta=0.2,     # Type II error rate (1 - power)
  L=0.05        # Half-width of the confidence interval
)

# Display the sample planning report
print(report)
```
For more detailed examples and usage, please refer to the help documentation.

## 🌐 Web Application
For users who prefer a graphical interface, the companion web app provides the same calculations with guided parameter entry and real‑time visualization.
Access the app here: https://spcal-demo.vercel.app/

