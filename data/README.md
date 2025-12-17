# Data Notes

This document describes the data environment, sources, variable definitions, and transformations used for the Bayesian estimation of the US–China DSGE model. It is intended to ensure transparency, reproducibility, and consistency with standard macroeconomic practice.

---

## 1. Target Environment

All data processing and estimation are conducted in the following environment:

- **Operating System:** Windows 11  
- **MATLAB Version:** R2023b  
- **Dynare Version:** 6.4  

All scripts are designed to be fully reproducible within this environment.

---

## 2. Data Sources (To Be Finalized)

Macroeconomic time series are obtained from standard international databases. Exact series identifiers, links, and vintages will be added once finalized.

### Federal Reserve Economic Data (FRED)
- US real GDP
- US CPI inflation
- Federal Funds Rate
- Additional US macroeconomic indicators as needed

### IMF – International Financial Statistics (IFS)
- China CPI
- China short-term interest rates
- Exchange rate data (USD/CNY)

### Bank for International Settlements (BIS)
- Financial indicators
- Cross-border banking and financial exposure measures
- Dollarization-related proxies (where applicable)

### World Bank
- National accounts data
- Trade and external sector indicators
- Supplementary macroeconomic series

---

## 3. Observable Variables

The following macroeconomic series are used as observables in the Bayesian estimation:

- Output (United States)
- Inflation (United States)
- Policy interest rate (United States)
- Output (China)
- Inflation (China)
- Policy interest rate (China)
- Exchange rate (USD/CNY)
- Dollarization proxy for China (where available)

Each observable is mapped directly to a corresponding theoretical variable in the DSGE model.

---

## 4. Series Identifiers and Units

Exact series identifiers, data links, and units of measurement will be documented here once finalized.

- **Series IDs:** To be added  
- **Frequency:** Quarterly  
- **Units:** As reported by source (levels, index values, or percentages)

---

## 5. Default Data Transformations

Unless otherwise stated, the following transformations are applied consistently across series.

### Inflation
- Computed as the log-difference of the price index  
- Annualized using a factor of 400:
  \[
  \pi_t = 400 \times \Delta \log(P_t)
  \]

### Output
- Measured as log real GDP  
- Used either as:
  - Output growth, or
  - Output gap (deviation from trend or steady state)  
- Exact specification is reported alongside estimation results

### Policy Interest Rate
- Short-term nominal interest rate  
- Expressed in quarterly frequency and annualized percentage terms

### Exchange Rate
- Nominal exchange rate (USD/CNY)  
- Used either in:
  - Levels, or
  - Log-differences  
- The chosen transformation is explicitly stated in the estimation section

### Trade Balance
- Defined as the ratio of net exports to GDP  
- Expressed in percentage terms

---

## 6. Data Treatment and Adjustments

- All series are seasonally adjusted when available.
- Missing observations are handled conservatively to preserve sample consistency.
- No ad hoc filtering is applied unless explicitly stated (e.g., HP filter for output gap construction).
- All transformations are implemented via scripted routines to ensure reproducibility.

---

## 7. Replication and Transparency

- All raw-to-processed data transformations are fully documented in code.
- Large raw datasets are not committed to the repository but can be downloaded from original sources.
- This documentation will be updated once final data choices are locked.

---

## Final Remark

These data notes provide a transparent foundation for the Bayesian estimation and IRF analysis. Any deviations from the default transformations or data definitions are explicitly documented in the estimation section of the paper.
