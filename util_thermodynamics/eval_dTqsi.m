function dTqsi = eval_dTqsi(T,p,qt)
%{
  Calculate partial derivative of (ice) saturation specific humidity with respect to temperature.

  Units: kg/kg/T.

  Args:
    T (double): air temperature in K.
    p (double): air pressure in Pa.
    qt (double): total water mixing ratio in kg/kg.
    EPS (constant): ratio of gas constant of dry air to gas constant of water vapor.

  Returns:
    double: derivative of saturation specific humidity w.r.t. T.
%}
  global EPS RV LS
  esi = eval_esi(T);
  % With es from Eq. (4.4.15) in Emanual (1994)
  %   accurate within <1e-3 rel. error for -30<=T<=0-deg-C
  % Rel. error increases exponentially with T
  if nargin==2 %assuming condensate-free
    dTqsi = EPS*LS/RV .* p .* esi ./ ( ( p-esi*(1-EPS) ).*T ).^2;
	else %given qt (suitable for plume computations)
    dTqsi = EPS*LS/RV .* (1-qt) .* p .* esi ./ ( ( p-esi ).*T ).^2;
  end