%{
Thermodynamic constants used for plume buoyancy calculations.

Notations (in uppercase) & values from Appendix 2 (p. 566-567) of Emanuel (1994).

Ref: Emanuel, K. A. 1994. Atmospheric Convection. Oxford University Press, 580 pp. 
%}
global RD RV CVD CPD CVV CPV CL CPVMCL EPS LV0 LS LF T0 g

% Units: J/Kg/K
RD = 287.04; % gas constant of dry air
RV = 461.5; % gas constant of water vapor
CVD = 719; % heat capacity at constant volume for dry air
CPD = 1005.7; % heat capacity at constant pressure for dry air
CVV = 1410; % heat capacity at constant volume of water vapor
CPV = 1870; % heat capacity at constant pressure of water vapor
CL = 4190; % heat capacity of liquid water (above freezing)
CPVMCL = CL-CPV; % CPVMCL seems to be a common notation for this value

% Units: dimensionless
EPS = RD/RV;

% Units: J/Kg
LV0 = 2.501e6; % latent heat of vaporization at 0-deg-C
LS = 2.834e6; % latent heat of sublimation (-100<=T<= 0-deg-C)
LF = 0.3337e6; % latent heat of fusion at 0-deg-C (lv0-ls)

% Units: K
% T0 = 273.16; % triple point
T0 = 273.1636783445389; % triple point so that eval_es(T0)=eval_esi(T0)

% Units: m/s^2
g = 9.80665; % standard gravity