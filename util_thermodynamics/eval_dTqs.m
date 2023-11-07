function dTqs = eval_dTqs(T,p,qt)
%{
  Calculate partial derivative of (liquid) saturation specific humidity with respect to temperature.

  Units: kg/kg/T.

  Args:
    T (double): air temperature in K.
    p (double): air pressure in Pa.
    qt (double): total water mixing ratio in kg/kg.
    EPS (constant): ratio of gas constant of dry air to gas constant of water vapor.

  Returns:
    double: derivative of saturation specific humidity w.r.t. T.
%}
  global EPS RV
  es = eval_es(T);
  % With es from Eq. (4.4.13) in Emanual (1994)
  %   accurate within <1e-3 rel. error for 0<=T<=20-deg-C
  %   <1.4e-3 up to T<=40-deg-C
  % Rel. error increases exponentially with T
  % dTqs = EPS.*p./( p-es*(1-EPS) ).^2 .* d(es)/dT
  if nargin==2 %assuming condensate-free
    dTqs = EPS/RV .* p .* lv(T) .* es ./ ( ( p-es*(1-EPS) ).*T ).^2;
 	else %given qt (suitable for plume computations)
    dTqs = EPS/RV .* (1-qt) .* p .* lv(T) .* es ./ ( ( p-es ).*T ).^2;
  end