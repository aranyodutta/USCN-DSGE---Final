// =====================================================
// CN-US Two-Country DSGE with Dollarization and UIP
// BAYESIAN ESTIMATION VERSION
// =====================================================

// -------------------------
// PARAMETERS
// -------------------------
parameters
  sigma eta beta theta_US theta_CN alpha varphi
  kappa_US kappa_CN omega_2 omega_4 psi varsigma delta
  sigma_omega oomega_2
  rho_US rho_CN rho_dol rho_mu
  phi_pi_US phi_pi_CN phi_x_US phi_x_CN phi_de_CN
  varrho_US varrho_CN
  delta_rUS gamma_dol
  chi_tb chi_dcp
  lambda_ss;

// -------------------------
// PARAMETER CALIBRATION
// -------------------------
sigma       = 4.5;
eta         = 2.5;
beta        = 0.99;
theta_US    = 0.9;
theta_CN    = 0.75;
alpha       = 0.6;
varphi      = 3;

kappa_US    = (1-theta_US)*(1-theta_US*beta)/theta_US;
kappa_CN    = (1-theta_CN)*(1-theta_CN*beta)/theta_CN;

omega_2     = alpha*2*(1-alpha)*(sigma*eta-1);
omega_4     = alpha*4*(1-alpha)*(sigma*eta-1);

psi         = (omega_4+1)*(1+varphi);
varsigma    = (omega_2+1)*sigma + (omega_4+1)*varphi;

delta =
  sigma^2*(2*omega_2+1)
  + 2*sigma*varphi*(omega_2+1)*(omega_4+1)
  + (omega_4+1)^2*varphi^2;

sigma_omega = (omega_2+1)*sigma/(omega_4+1);
oomega_2    = alpha*sigma*(1+varphi)*(varsigma + omega_2*sigma)/delta;

rho_US     = 0.6;
rho_CN     = 0.6;
rho_dol    = 0.7;
rho_mu     = 0.7;

phi_pi_US  = 1.5;
phi_pi_CN  = 1.5;
phi_x_US   = 0.5;
phi_x_CN   = 0.5;
phi_de_CN  = 0.3;

varrho_US  = 0.4;
varrho_CN  = 0.4;
delta_rUS  = 0.35;
gamma_dol  = 0.5;

chi_tb     = 1.5;
chi_dcp    = 0.4;
lambda_ss  = 0.3;

// -------------------------
// ENDOGENOUS VARIABLES
// -------------------------
var
  x_US pi_US r_US pi_US_CPI
  x_CN pi_CN r_CN pi_CN_CPI
  rbar_US rbar_CN
  a_US a_CN
  mc_US mc_CN
  y_US y_CN ybar_US ybar_CN
  p_US p_CN e
  dol_CN mu
  tb_US tb_CN
  gdp_US gdp_CN
  tbgdp_US tbgdp_CN
  de
  // Observed variables (with measurement error)
  dy dinfl dr dy_c dinfl_c dr_c;

// -------------------------
// EXOGENOUS SHOCKS
// -------------------------
varexo
  m_US m_CN
  xi_US xi_CN
  xi_dol
  e1 e2 e3 e4 e5 e6;

// -------------------------
// MODEL EQUATIONS
// -------------------------
model(linear);

  // IS CURVES
  x_US =
    x_US(+1)
    - (omega_4+1)/((omega_2+1)*sigma)*r_US
    + (omega_4+1)/((omega_2+1)*sigma)*pi_US(+1)
    + omega_2/(omega_2+1)*(x_CN(+1)-x_CN)
    + (omega_4+1)/((omega_2+1)*sigma)*rbar_US;

  x_CN =
    x_CN(+1)
    - (omega_4+1)/((omega_2+1)*sigma)*r_CN
    + (omega_4+1)/((omega_2+1)*sigma)*pi_CN(+1)
    + omega_2/(omega_2+1)*(x_US(+1)-x_US)
    + (omega_4+1)/((omega_2+1)*sigma)*rbar_CN;

  // PHILLIPS CURVES
  pi_US =
    beta*pi_US(+1)
    + kappa_US*varsigma/(omega_4+1)*x_US
    + kappa_US*omega_2*sigma/(omega_4+1)*x_CN
    + kappa_US*r_US;

  pi_CN =
    beta*pi_CN(+1)
    + kappa_CN*varsigma/(omega_4+1)*x_CN
    + kappa_CN*omega_2*sigma/(omega_4+1)*x_US
    + kappa_CN*r_CN;

  // MONETARY POLICY RULES
  r_US =
    varrho_US*r_US(-1)
    + (1-varrho_US)*(phi_pi_US*pi_US_CPI + phi_x_US*x_US)
    + m_US;

  r_CN =
    varrho_CN*r_CN(-1)
    + (1-varrho_CN)*(phi_pi_CN*pi_CN_CPI + phi_x_CN*x_CN + phi_de_CN*de + delta_rUS*r_US)
    + m_CN;

  // BACKWARD-LOOKING UIP
  r_CN - r_US = de + mu;

  // DOLLARIZATION AND RISK PREMIUM
  mu     = rho_mu*mu(-1) + gamma_dol*dol_CN;
  dol_CN = rho_dol*dol_CN(-1) + xi_dol;

  // CPI INFLATION
  pi_US_CPI =
    pi_US
    + alpha*sigma/(omega_4+1)*(x_US - x_US(-1) - x_CN + x_CN(-1))
    + oomega_2*(a_US - a_US(-1) - a_CN + a_CN(-1));

  pi_CN_CPI =
    pi_CN
    + alpha*sigma/(omega_4+1)*(x_CN - x_CN(-1) - x_US + x_US(-1))
    + oomega_2*(a_CN - a_CN(-1) - a_US + a_US(-1))
    + chi_dcp*lambda_ss*de;

  // MARGINAL COSTS
  mc_US = varsigma/(omega_4+1)*x_US + omega_2*sigma/(omega_4+1)*x_CN + r_US;
  mc_CN = varsigma/(omega_4+1)*x_CN + omega_2*sigma/(omega_4+1)*x_US + r_CN;

  // OUTPUT AND NATURAL RATES
  ybar_US = varsigma*psi/delta*a_US - omega_2*sigma*psi/delta*a_CN;
  ybar_CN = varsigma*psi/delta*a_CN - omega_2*sigma*psi/delta*a_US;

  y_US = x_US + ybar_US;
  y_CN = x_CN + ybar_CN;

  rbar_US =
    - sigma*(1-rho_US)*psi*((omega_2+1)*varsigma - omega_2^2*sigma)/((omega_4+1)*delta)*a_US
    - sigma*(1-rho_US)*omega_2*psi*(varsigma - sigma*(omega_2+1))/((omega_4+1)*delta)*a_CN;

  rbar_CN =
    - sigma*(1-rho_CN)*psi*((omega_2+1)*varsigma - omega_2^2*sigma)/((omega_4+1)*delta)*a_CN
    - sigma*(1-rho_CN)*omega_2*psi*(varsigma - sigma*(omega_2+1))/((omega_4+1)*delta)*a_US;

  // TRADE BALANCE
  tb_US = -chi_tb*e + alpha*(y_CN - y_US);
  tb_CN =  chi_tb*e + alpha*(y_US - y_CN);

  gdp_US = y_US + p_US;
  gdp_CN = y_CN + p_CN;

  tbgdp_US = tb_US - gdp_US;
  tbgdp_CN = tb_CN - gdp_CN;

  // PRICE LEVELS
  p_US = p_US(-1) + pi_US_CPI;
  p_CN = p_CN(-1) + pi_CN_CPI;

  // PRODUCTIVITY SHOCKS
  a_US = rho_US*a_US(-1) + xi_US;
  a_CN = rho_CN*a_CN(-1) + xi_CN;

  // EXCHANGE RATE GROWTH
  de = e - e(-1);

  // -------------------------
  // MEASUREMENT EQUATIONS
  // -------------------------
  dy      = y_US - y_US(-1) + e1;
  dinfl   = p_US - p_US(-1) + e2;
  dr      = r_US - r_US(-1) + e3;
  dy_c    = y_CN - y_CN(-1) + e4;
  dinfl_c = p_CN - p_CN(-1) + e5;
  dr_c    = r_CN - r_CN(-1) + e6;

end;

// -------------------------
// INITIAL VALUES
// -------------------------
initval;
  x_US = 0; pi_US = 0; r_US = 0; pi_US_CPI = 0;
  x_CN = 0; pi_CN = 0; r_CN = 0; pi_CN_CPI = 0;
  rbar_US = 0; rbar_CN = 0;
  a_US = 0; a_CN = 0;
  mc_US = 0; mc_CN = 0;
  y_US = 0; y_CN = 0; ybar_US = 0; ybar_CN = 0;
  p_US = 0; p_CN = 0; e = 0; de = 0;
  dol_CN = 0; mu = 0;
  tb_US = 0; tb_CN = 0;
  gdp_US = 0; gdp_CN = 0;
  tbgdp_US = 0; tbgdp_CN = 0;
  dy = 0; dinfl = 0; dr = 0;
  dy_c = 0; dinfl_c = 0; dr_c = 0;
end;

steady;
check;

// -------------------------
// SHOCKS
// -------------------------
shocks;
  var m_US;      stderr 0.25;
  var m_CN;      stderr 0.25;
  var xi_US;     stderr 1.0;
  var xi_CN;     stderr 1.0;
  var xi_dol;    stderr 0.1;
  var e1;        stderr 0.1;
  var e2;        stderr 0.1;
  var e3;        stderr 0.1;
  var e4;        stderr 0.1;
  var e5;        stderr 0.1;
  var e6;        stderr 0.1;
end;

// -------------------------
// OBSERVABLES
// -------------------------
varobs dy dinfl dr dy_c dinfl_c dr_c;

// =====================================================
// ESTIMATED PARAMETERS
// =====================================================
estimated_params;
  // AR(1) persistence
  rho_US,     beta_pdf,       0.6,  0.2;
  rho_CN,     beta_pdf,       0.6,  0.2;
  rho_dol,    beta_pdf,       0.7,  0.15;
  rho_mu,     beta_pdf,       0.7,  0.15;
  
  // Price stickiness
  theta_US,   beta_pdf,       0.75, 0.1;
  theta_CN,   beta_pdf,       0.75, 0.1;
  
  // Risk aversion
  sigma,      gamma_pdf,      4.5,  0.5;
  
  // Taylor rule
  phi_pi_US,  normal_pdf,     1.5,  0.3;
  phi_pi_CN,  normal_pdf,     1.5,  0.3;
  phi_x_US,   normal_pdf,     0.5,  0.2;
  phi_x_CN,   normal_pdf,     0.5,  0.2;
  phi_de_CN,  normal_pdf,     0.3,  0.1;
  
  // Interest rate smoothing
  varrho_US,  beta_pdf,       0.4,  0.2;
  varrho_CN,  beta_pdf,       0.4,  0.2;
  
  // Spillovers
  delta_rUS,  normal_pdf,     0.35, 0.1;
  gamma_dol,  normal_pdf,     0.5,  0.2;
  
  // DCP and trade
  lambda_ss,  beta_pdf,       0.3,  0.15;
  chi_tb,     gamma_pdf,      1.5,  0.3;
  
  // Structural shocks
  stderr m_US,    inv_gamma_pdf, 0.25, inf;
  stderr m_CN,    inv_gamma_pdf, 0.25, inf;
  stderr xi_US,   inv_gamma_pdf, 1.0,  inf;
  stderr xi_CN,   inv_gamma_pdf, 1.0,  inf;
  stderr xi_dol,  inv_gamma_pdf, 0.1,  inf;
  
  // Measurement errors
  stderr e1,      inv_gamma_pdf, 0.1,  inf;
  stderr e2,      inv_gamma_pdf, 0.1,  inf;
  stderr e3,      inv_gamma_pdf, 0.1,  inf;
  stderr e4,      inv_gamma_pdf, 0.1,  inf;
  stderr e5,      inv_gamma_pdf, 0.1,  inf;
  stderr e6,      inv_gamma_pdf, 0.1,  inf;
end;

// =====================================================
// BAYESIAN ESTIMATION
// =====================================================
estimation(
  datafile = data_6,
  first_obs = 1,
  nobs = 93,
  mh_replic = 20000,
  mh_nblocks = 2,
  mh_drop = 0.3,
  mh_jscale = 0.2,
  diffuse_filter,
  mode_check,
  mode_compute = 6,
  bayesian_irf,
  plot_priors = 1
);