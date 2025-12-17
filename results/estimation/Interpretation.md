# How to Read and Interpret the Bayesian Estimation

This section explains the **logic, workflow, and interpretation strategy** used in the Bayesian estimation of the two-country DSGE model. The objective is to clarify how priors, data, and the model jointly determine the posterior estimates and how these estimates should be interpreted in the context of the economic mechanisms studied in the paper.

---

## 1. Purpose of the Bayesian Estimation

Bayesian estimation is used to:
- Discipline the DSGE model with observed macroeconomic data
- Quantify uncertainty around structural parameters
- Evaluate the empirical plausibility of the model’s transmission mechanisms
- Provide a foundation for impulse response analysis

The estimation combines **prior economic beliefs** with **information from the data** to produce posterior distributions for model parameters.

---

## 2. Step 1: Selection of Observable Variables

The estimation begins with the choice of observable macroeconomic time series. These variables are selected to correspond closely to theoretical objects in the model.

Typical observables include:
- Output growth
- Inflation rates
- Short-term interest rates
- Dollarization proxies for China

Each observable is mapped directly to a model variable, ensuring consistency between theory and data.

---

## 3. Step 2: Specification of Priors

Priors encode **economic knowledge and empirical regularities** documented in the literature.

### Key Principles
- Priors are centered on standard values used in New Keynesian DSGE models
- Prior distributions reflect parameter domains (e.g., beta, normal, gamma)
- Prior variances allow sufficient flexibility for the data to be informative

Examples:
- High persistence priors for policy rates and dollarization
- Taylor-rule coefficients centered above unity for inflation
- Moderate price stickiness consistent with quarterly data

Priors are not intended to dominate the results but to regularize estimation.

---

## 4. Step 3: Likelihood and State-Space Representation

The log-linearized DSGE model is written in **state-space form**, allowing the use of the Kalman filter.

- The **measurement equation** links observed data to model variables
- The **transition equation** captures the model’s dynamic structure

The likelihood measures how well the model explains the observed data given a set of parameters.

---

## 5. Step 4: Posterior Distribution

The posterior distribution combines:
- The prior distribution
- The likelihood function

Bayes’ rule updates beliefs about parameters based on the information contained in the data.

Posterior results are summarized using:
- Posterior means
- Highest posterior density (HPD) intervals
- Marginal posterior distributions

These summaries form the basis for economic interpretation.

---

## 6. Step 5: Interpretation of Structural Parameters

Posterior estimates are interpreted in light of economic theory:

- **Preference parameters** inform consumption and labor responses
- **Price stickiness parameters** determine inflation persistence
- **Policy rule coefficients** measure central bank aggressiveness
- **Dollarization parameters** quantify how financial dollarization alters monetary transmission

Consistency with the literature strengthens confidence in the model.

---

## 7. Step 6: Shock Processes and Volatilities

Estimated shock processes provide insight into:
- The persistence of structural shocks
- The relative importance of monetary vs. real disturbances
- The magnitude of dollarization shocks in China

High estimated persistence implies prolonged macroeconomic effects, which directly influence IRFs.

---

## 8. Step 7: Model Fit and Data Density

Model fit is evaluated using:
- The log data density
- Comparison of prior vs. posterior distributions
- Visual inspection of filtered and smoothed variables

A higher data density indicates better empirical performance relative to alternative specifications.

---

## 9. Step 8: MCMC Diagnostics and Convergence

Reliable inference requires well-behaved Markov Chain Monte Carlo (MCMC) sampling.

Diagnostics include:
- Trace plots
- Autocorrelation functions
- Inefficiency factors
- Brooks–Gelman convergence statistics

Convergence ensures that posterior summaries reflect the true posterior distribution rather than sampling artifacts.

---

## 10. Step 9: Link Between Estimation and IRFs

Bayesian estimation directly shapes the impulse response analysis:

- IRFs are computed using posterior mean parameters
- Uncertainty in parameters translates into uncertainty in IRFs
- Estimated policy rules and shock processes determine IRF magnitudes and persistence

Thus, IRFs should always be read as **model-consistent and data-informed responses**.

---

## 11. Interpretation Strategy

To interpret the Bayesian estimation correctly:
1. Start with the **priors** and their economic motivation
2. Examine **posterior means and HPD intervals**
3. Compare **prior vs. posterior shifts**
4. Assess **parameter plausibility and identification**
5. Link estimates to **IRF behavior**

This approach ensures disciplined and transparent inference.

---

## Final Remark

The Bayesian estimation is not an end in itself, but a tool that grounds the DSGE model in data. It ensures that the impulse response analysis reflects empirically plausible dynamics while preserving the structural interpretation necessary for policy analysis.
