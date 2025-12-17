# US–China DSGE Model  
**A Two-Country New Keynesian DSGE Framework with Dollarization**

This repository contains a two-country Dynamic Stochastic General Equilibrium (DSGE) model of the United States and China designed to study **dollar dominance**, **financial dollarization**, and **de-dollarization dynamics**. The model is implemented in **Dynare/MATLAB** and produces impulse response functions (IRFs), simulated moments, and supporting documentation for academic research.

---

## 1. Research Motivation

The dominance of the US dollar in global financial markets shapes monetary policy transmission, exchange rate dynamics, and macroeconomic spillovers. Emerging economies with partial financial dollarization—such as China—face constraints on monetary autonomy that alter the effectiveness of domestic policy.

This project develops a structural two-country New Keynesian DSGE model to analyze:
- How US monetary policy shocks transmit internationally
- How financial dollarization affects monetary policy effectiveness in China
- The macroeconomic consequences of de-dollarization
- Cross-border spillovers in inflation, output, and interest rates

---

## 2. Model Overview

### Economic Structure
- **Home economy:** United States (advanced economy, issuer of the global reserve currency)
- **Foreign economy:** China (large emerging economy with partial financial dollarization)

Each country features:
- Representative households with CES consumption baskets
- Monopolistically competitive firms with Calvo price stickiness
- New Keynesian Phillips curves
- Taylor-type monetary policy rules

### Key Features
- Open-economy trade linkages and terms-of-trade effects
- Nominal exchange rate dynamics
- Explicit dollarization channel in China’s monetary policy rule
- Spillover effects from US monetary policy to China

---

## 3. Shocks Analyzed

The model focuses on two primary sources of uncertainty:

1. **US Monetary Policy Shock**
   - Contractionary or expansionary shocks to the Federal Reserve’s policy rate
   - Transmission through interest rates, inflation, output, and the exchange rate

2. **China Dollarization / De-dollarization Shock**
   - Shocks to the degree of financial dollarization in China
   - Effects on monetary policy transmission, inflation dynamics, and output
   - Spillovers to the US economy

Impulse response functions (IRFs) are the main tool for analyzing these dynamics.

---

## 4. Software Requirements

To run the model and reproduce the results, you need:

- **MATLAB**  
  - Tested on R2023b (earlier versions may work)
- **Dynare 5.x**  
  - Must be installed and added to the MATLAB path

---

## 5. Quick Start Guide

To generate all baseline results and figures:

1. Open MATLAB  
2. Navigate to the project root directory  
3. Run:
   ```matlab
   cd code/scripts
   run_all.m
