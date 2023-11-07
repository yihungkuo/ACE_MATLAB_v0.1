function qs = eval_qs(T, p, qt)
%{
  Calculate saturation specific humidity with respect to liquid water.

  Units: kg/kg (unitless).

  Args:
    T (double): air temperature in K.
    p (double): air pressure in Pa.
    qt (double): total water mixing ratio in kg/kg.
    EPS (constant): ratio of gas constant of dry air to gas constant of water vapor.

  Returns:
    double: saturation specific humidity.
%}
%   global EPS
  es = eval_es(T);
  
  if nargin==2 %assuming condensate-free
    qs = eval_q(es,p);
  else %given qt (suitable for plume computations)
    qs = eval_q(es,p,qt);
  end