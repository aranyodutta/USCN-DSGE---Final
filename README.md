# US–China DSGE Model
## A Two-Country New Keynesian DSGE Framework with Dollarization

This repository contains a two-country Dynamic Stochastic General Equilibrium (DSGE) model (United States and China) implemented in Dynare/MATLAB. It is designed to study international monetary spillovers and the macroeconomic implications of financial dollarization and de-dollarization. The code produces impulse response functions (IRFs), simulated moments, and supporting documentation.

**Literature note.** The model’s open-economy two-country New Keynesian core follows standard frameworks in the literature and is **adapted in part from Jang and Okano (2013)**, with additional extensions to study dollarization and policy spillovers (see **References**). 

---

## 1. Research Motivation
US dollar dominance in global finance can shape monetary policy transmission, exchange rate dynamics, and cross-border macroeconomic spillovers. In economies with meaningful foreign-currency exposure, partial financial dollarization can constrain monetary autonomy and alter domestic policy effectiveness.

This project develops and extends a two-country New Keynesian DSGE framework to analyze:

- How US monetary policy shocks transmit internationally
- How financial dollarization can affect monetary policy effectiveness in China
- The macroeconomic consequences of de-dollarization
- Cross-border spillovers in inflation, output, and interest rates

---

## 2. Model Overview

### Economic Structure
- **Home economy:** United States (advanced economy; issuer of the dominant reserve currency)
- **Foreign economy:** China (large open economy with partial financial dollarization)

Each country features:
- Representative households with CES consumption baskets
- Monopolistically competitive firms with Calvo price stickiness
- New Keynesian Phillips curves
- Taylor-type monetary policy rules

### Key Features
- Open-economy trade linkages and terms-of-trade effects
- Nominal exchange rate dynamics
- A dollarization state in China that enters the monetary policy rule (dollarization / de-dollarization channel)
- Spillover effects from US monetary policy to China

---

## 3. Shocks Analyzed
The model focuses on two primary sources of uncertainty:

1) **US Monetary Policy Shock**
- Contractionary or expansionary shocks to the Federal Reserve’s policy rate
- Transmission through interest rates, inflation, output, and the exchange rate

2) **China Dollarization / De-dollarization Shock**
- Shocks to the degree of financial dollarization in China (modeled as a disturbance to the dollarization state)
- Effects on monetary policy transmission, inflation dynamics, and output
- Spillovers back to the US economy

Impulse response functions (IRFs) are the main tool for analyzing these dynamics.

---

## 4. Software Requirements
To run the model and reproduce results, you need:

- **MATLAB** (tested on R2023b; earlier versions may work)
- **Dynare 5.x** (installed and added to the MATLAB path)

---

## 5. Quick Start Guide
To generate baseline results and figures:

1. Open MATLAB  
2. Navigate to the project root directory  
3. Run:

```matlab
cd code/scripts
run_all.m
