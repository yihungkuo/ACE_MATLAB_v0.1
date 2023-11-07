function q = eval_q(e, p, qt)
%{
  Calculate specific humidity.

  Units: kg/kg (unitless).

  Args:
    e (double): vapor pressure in Pa.
    p (double): air pressure in Pa.
    qt (double): total water mixing ratio in kg/kg.
    EPS (constant): ratio of gas constant of dry air to gas constant of water vapor.

  Returns:
    double: specific humidity.
%}  
  global EPS
  if nargin==2 %assuming condensate-free
    q = EPS * e ./ ( p-e*(1-EPS) );
  else %given qt (suitable for plume computations)
    q = EPS * e .* (1-qt) ./ (p-e);
  end