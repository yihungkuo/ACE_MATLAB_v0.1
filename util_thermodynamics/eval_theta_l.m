function theta_l = eval_theta_l(T, q, ql, p, P0)
%{
  Calculate liquid water potential temperature.

  Units: K.

  Args:
    T (double): air temperature in K.
    q (double): specific humidity of water vapor in kg/kg.
    ql (double): specific liquid water content in kg/kg.
    p (double): air pressure in Pa.
    RD (constant): gas constant of dry air in J/Kg/K.
    RV (constant): gas constant of water vapor in J/Kg/K.
    CPD (constant): heat capacity at constant pressure for dry air in J/Kg/K.
    CPV (constant): heat capacity at constant pressure for water vapor in J/Kg/K.
    EPS (constant): ratio of gas constant of dry air to gas constant of water vapor.
    P0 (constant): reference pressure (1,000 hPa) in Pa.

  Returns:
    double: liquid water potential temperature following Eq. (2.44) in Stevens & Siebesma (2020)
              which is equivalent to Eq. (4.5.15) in Emanuel (1994)
%}
  global RD RV CPD CPV EPS
  if nargin==4
    P0 = 1e5;
  end
  LV = lv(T);
  
  qt = q + ql;
  
  RL = RD*( 1 + (1-EPS)/EPS*qt );
  R = (1-qt)*RD + q*RV;
  
  CPL = CPD + qt*(CPV-CPD);
  
  chi = RL./CPL;
  gamma = RV*qt ./ CPL;
  
  theta_l = T .* (P0*R./p./RL).^chi ...
              .* (qt./q).^gamma ...
              .* exp( -LV.*ql./CPL./T );

%% Below: Eq. (4.5.15) from Emanuel (1994); equivalent to Eq. (2.44) in Stevens & Siebesma (2020)
%   r = q ./ (1-q-ql);
%   rl = ql ./ (1-q-ql);
%   rt = r+rl;
% 
%   CPL = CPD + rt*CPV;
% 
%   chi = (RD + rt*RV) ./ CPL;
%   gamma = rt*RV ./ CPL;
%   
%   theta_l = T .* (P0./p).^chi ...
%               .* ( 1-rl./(EPS+rt) ).^chi ...
%               .* (1-rl./rt).^(-gamma) ...
%               .* exp(-LV.*rl./CPL./T);
            
  