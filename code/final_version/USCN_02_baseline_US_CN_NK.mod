var x_US pi_US r_US pi_US_CPI 
    x_CN pi_CN r_CN pi_CN_CPI 
    rbar_US rbar_CN 
    a_US a_CN 
    mc_US mc_CN 
    y_US y_CN ybar_US ybar_CN 
    p_US p_CN e s q;

varexo m_US m_CN xi_US xi_CN;

parameters sigma eta beta theta_US theta_CN alpha varphi 
    kappa_US kappa_CN omega_2 omega_4 psi varsigma delta 
    sigma_omega oomega_2 rho_US rho_CN 
    phi_pi_US phi_pi_CN phi_x_US phi_x_CN varrho_US varrho_CN;

sigma = 4.5;
eta   = 2.5;
beta  = 0.99;

theta_US = 0.9;
theta_CN = 0.75;
alpha    = 0.6;
varphi   = 3;

kappa_US = (1 - theta_US) * (1 - theta_US * beta) / theta_US;
kappa_CN = (1 - theta_CN) * (1 - theta_CN * beta) / theta_CN;

omega_2 = alpha * 2 * (1 - alpha) * (sigma * eta - 1);
omega_4 = alpha * 4 * (1 - alpha) * (sigma * eta - 1);
psi       = (omega_4 + 1) * (1 + varphi);
varsigma  = (omega_2 + 1) * sigma + (omega_4 + 1) * varphi;
delta     = sigma^2 * (2 * omega_2 + 1) + 2 * sigma * varphi * (omega_2 + 1) * (omega_4 + 1) + (omega_4 + 1)^2 * varphi^2;
sigma_omega = (omega_2 + 1) * sigma / (omega_4 + 1);
oomega_2    = alpha * sigma * (1 + varphi) * (varsigma + omega_2 * sigma) / delta;

rho_US = 0;      
rho_CN = 0.55;

phi_pi_US      = 1.5;
phi_pi_CN      = 1.5;
phi_x_US       = 0.5;
phi_x_CN       = 0.5;
varrho_US      = 0.4;
varrho_CN      = 0.4;

model(linear);

  // -------------------------
  // US IS curve
  // -------------------------
  x_US = x_US(+1) 
         - (omega_4 + 1)/((omega_2 + 1)*sigma) * r_US 
         + (omega_4 + 1)/((omega_2 + 1)*sigma) * pi_US(+1)
         + omega_2/(omega_2 + 1) * x_CN(+1) - omega_2/(omega_2 + 1) * x_CN
         + (omega_4 + 1)/((omega_2 + 1)*sigma) * rbar_US;

  // -------------------------
  // US Phillips curve
  // -------------------------
  pi_US = beta * pi_US(+1) 
          + kappa_US * varsigma/(omega_4 + 1) * x_US
          + kappa_US * omega_2 * sigma/(omega_4 + 1) * x_CN 
          + kappa_US * r_US;

  // -------------------------
  // China IS curve
  // -------------------------
  x_CN = x_CN(+1) 
         - (omega_4 + 1)/((omega_2 + 1)*sigma) * r_CN
         + (omega_4 + 1)/((omega_2 + 1)*sigma) * pi_CN(+1)
         + omega_2/(omega_2 + 1) * x_US(+1) - omega_2/(omega_2 + 1) * x_US
         + (omega_4 + 1)/((omega_2 + 1)*sigma) * rbar_CN;

  // -------------------------
  // China Phillips curve
  // -------------------------
  pi_CN = beta * pi_CN(+1) 
          + kappa_CN * varsigma/(omega_4 + 1) * x_CN
          + kappa_CN * omega_2 * sigma/(omega_4 + 1) * x_US 
          + kappa_CN * r_CN;

  // -------------------------
  // Monetary policy rules (Taylor rules)
  // -------------------------
  r_US = varrho_US * r_US(-1) + (1 - varrho_US) * (phi_pi_US * pi_US_CPI + phi_x_US * x_US) + m_US;
  r_CN = varrho_CN * r_CN(-1) + (1 - varrho_CN) * (phi_pi_CN * pi_CN_CPI + phi_x_CN * x_CN) + m_CN;

  // -------------------------
  // Natural interest rates
  // -------------------------
  rbar_US = - sigma * (1 - rho_US) * psi * ((omega_2 + 1) * varsigma - omega_2^2 * sigma) / ((omega_4 + 1) * delta) * a_US
            - sigma * (1 - rho_US) * omega_2 * psi * (varsigma - sigma * (omega_2 + 1)) / ((omega_4 + 1) * delta) * a_CN;

  rbar_CN = - sigma * (1 - rho_US) * psi * ((omega_2 + 1) * varsigma - omega_2^2 * sigma) / ((omega_4 + 1) * delta) * a_CN
            - sigma * (1 - rho_US) * omega_2 * psi * (varsigma - sigma * (omega_2 + 1)) / ((omega_4 + 1) * delta) * a_US;

  // -------------------------
  // CPI inflation definitions
  // -------------------------
  pi_US_CPI = pi_US
              + alpha * sigma/(omega_4 + 1) * x_US - alpha * sigma/(omega_4 + 1) * x_US(-1)
              - alpha * sigma/(omega_4 + 1) * x_CN + alpha * sigma/(omega_4 + 1) * x_CN(-1)
              + oomega_2 * a_US - oomega_2 * a_US(-1) - oomega_2 * a_CN + oomega_2 * a_CN(-1);

  pi_CN_CPI = pi_CN
              + alpha * sigma/(omega_4 + 1) * x_CN - alpha * sigma/(omega_4 + 1) * x_CN(-1)
              - alpha * sigma/(omega_4 + 1) * x_US + alpha * sigma/(omega_4 + 1) * x_US(-1)
              + oomega_2 * a_CN - oomega_2 * a_CN(-1) - oomega_2 * a_US + oomega_2 * a_US(-1);

  // -------------------------
  // Marginal costs
  // -------------------------
  mc_US = varsigma/(omega_4 + 1) * x_US + omega_2 * sigma/(omega_4 + 1) * x_CN + r_US;
  mc_CN = varsigma/(omega_4 + 1) * x_CN + omega_2 * sigma/(omega_4 + 1) * x_US + r_CN;

  // -------------------------
  // Natural output
  // -------------------------
  ybar_US =  varsigma * psi / delta * a_US - omega_2 * sigma * psi / delta * a_CN;
  ybar_CN =  varsigma * psi / delta * a_CN - omega_2 * sigma * psi / delta * a_US;

  // -------------------------
  // Output
  // -------------------------
  y_US = x_US + ybar_US;
  y_CN = x_CN + ybar_CN;

  // -------------------------
  // Price levels
  // -------------------------
  p_US = pi_US_CPI + p_US(-1);
  p_CN = pi_CN_CPI + p_CN(-1);

  // -------------------------
  // Terms of trade, exchange rate
  // -------------------------
  s = sigma/(omega_4 + 1) * y_US - sigma/(omega_4 + 1) * y_CN;
  q = (1 - 2 * alpha) * s;
  e = q - p_CN + p_US;

  // -------------------------
  // Productivity shocks
  // -------------------------
  a_US = rho_US * a_US(-1) + xi_US;
  a_CN = rho_CN * a_CN(-1) + xi_CN;

end;

initval;
  x_US = 0; pi_US = 0; r_US = 0; pi_US_CPI = 0;
  x_CN = 0; pi_CN = 0; r_CN = 0; pi_CN_CPI = 0;
  rbar_US = 0; rbar_CN = 0;
  a_US = 0; a_CN = 0; xi_US = 0; xi_CN = 0;
  mc_US = 0; mc_CN = 0;
  y_US = 0; y_CN = 0; ybar_US = 0; ybar_CN = 0;
  p_US = 0; p_CN = 0; s = 0; q = 0; e = 0;
end;

steady;
check;

shocks;
  var xi_CN; stderr 1;   // China productivity shock
end;

stoch_simul(order=1, irf=40);
