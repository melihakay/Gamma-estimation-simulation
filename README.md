# Gamma-estimation-simulation

# Comparison of MoM and MLE for Gamma Distribution Parameters

## 📌 Project Overview
This project compares the performance of two statistical estimation methods, **Method of Moments (MoM)** and **Maximum Likelihood Estimation (MLE)**, for estimating the shape ($k$) and scale ($\theta$) parameters of the Gamma distribution.

The study involves theoretical derivations, application on a simulated dataset, and a comprehensive **Monte Carlo Simulation study** to evaluate finite-sample properties such as Bias, Variance, and Mean Squared Error (MSE).

## 📂 Repository Structure
* **`term_project.Rmd`**: The main RMarkdown source file containing all code, derivation (LaTeX), and analysis.
* **`term_project.pdf`**: The final compiled report presenting the results and visualizations.
* **`project_codes.R`**: A standalone R script containing the reproducible code for data generation and simulation (extracted from the report).
* **`README.md`**: Project documentation.

## ⚙️ Methodology
The project follows these main steps:
1.  **Theoretical Derivation:** Deriving the equations for MoM and MLE for Gamma parameters.
2.  **Point Estimation:** Applying the methods to a single generated dataset ($n=100, k=3, \theta=2$) to visualize fit using Histograms.
3.  **Monte Carlo Simulation:**
    * **Replications:** $R = 2000$
    * **Sample Sizes:** $n \in \{20, 50, 100\}$
    * **Scenarios:**
        * Scenario 1 (High Skewness): $k=1, \theta=2$
        * Scenario 2 (Moderate Skewness): $k=5, \theta=1$

## 📊 Key Findings
Based on the simulation results:
* **Consistency:** Both estimators are consistent; Bias and MSE decrease as sample size ($n$) increases.
* **Efficiency:** **MLE outperforms MoM** in terms of MSE in almost all scenarios.
* **Skewness Effect:** The advantage of MLE is most pronounced in highly skewed distributions (Scenario 1) and smaller sample sizes.

## 🚀 How to Run
To reproduce the results, you need **R** and **RStudio**.

1.  Clone this repository.
2.  Open `term_project.Rmd` in RStudio.
3.  Install required packages if missing:
    ```r
    install.packages(c("knitr", "rmarkdown"))
    ```
4.  Click the **Knit** button to generate the PDF report.

## 🛠 Technologies Used
* **R Language** (Statistical Computing)
* **RMarkdown** (Reproducible Reporting)
* **LaTeX** (Mathematical Typesetting)

## 👤 Author
[Rubar Akyıldız]
[akyildiz.rubar@metu.edu.tr]
