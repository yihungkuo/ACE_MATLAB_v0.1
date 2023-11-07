function theta_il = eval_theta_il(T, q, ql, qi, p, P0)
%{ 
  Units: K.

  This is the CONSERVED variable in the entraining plume calculation.

  Following Eq. (2.44) in Stevens & Siebesma (2020)
    which, with ice correction, is equivalent to Eq. (4.5.15) in Emanuel (1994)

  Args:
    T (double): air temperature in K.
    q (double): specific humidity of water vapor in kg/kg.
    ql (double): specific liquid water content in kg/kg.
    qi (double): specific ice content in kg/kg.
    p (double): air pressure in Pa.
    RD (constant): gas constant of dry air in J/Kg/K.
    RV (constant): gas constant of water vapor in J/Kg/K.
    CPD (constant): heat capacity at constant pressure for dry air in J/Kg/K.
    CPV (constant): heat capacity at constant pressure for water vapor in J/Kg/K.
    LS (constant): latent heat of sublimation in J/Kg/K (-100<=T<= 0-deg-C)
    EPS (constant): ratio of gas constant of dry air to gas constant of water vapor.
    P0 (constant): reference pressure (1,000 hPa) in Pa.

  Returns:
    double: ice-liquid water equivalent potential temperature.
%}
  global RD RV CPD CPV LS EPS
  if nargin==5
    P0 = 1e5;
  end
  LV = lv(T);
  
  qt = q + ql + qi;
  
  RL = RD*( 1 + (1-EPS)/EPS*qt );
  R = (1-qt)*RD + q*RV;
  
  CPL = CPD + qt*(CPV-CPD);
  
  chi = RL./CPL;
  gamma = RV*qt ./ CPL;
  
  theta_il = T .* (P0*R./p./RL).^chi ...
               .* (qt./q).^gamma ...
               .* exp( -( LV.*ql + LS*qi  )./CPL./T );
             
%% Below: Eq. (4.5.15) from Emanuel (1994); equivalent to Eq. (2.44) in Stevens & Siebesma (2020)  
%   r = q ./ (1-q-ql-qi);
%   rl = ql ./ (1-q-ql-qi);
%   ri = qi ./ (1-q-ql-qi);
%   rt = r+rl+ri;
% 
%   CPL = CPD + rt*CPV;
% 
%   chi = (RD + rt*RV) ./ CPL;
%   gamma = rt*RV ./ CPL;
% 	
%   theta_il = T .* (P0./p).^chi ...
%                .* ( 1-(rl+ri)./(EPS+rt) ).^chi ...
%                .* ( 1-(rl+ri)./rt ).^(-gamma) ...
%                .* exp((-LV.*rl-LS*ri)./CPL./T); 