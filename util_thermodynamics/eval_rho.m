function rho = eval_rho(T, q, ql, p)
%{
  Calculate air density.

  Units: Kg/m^3.

  Virtual (density) temperature is used in the calculation.

  Args:
    T (double): air temperature in K.
    q (double): specific humidity of water vapor in kg/kg.
    ql (double): specific liquid water content in kg/kg.
    p (double): air pressure in Pa.
    RD (constant): gas constant of dry air in J/Kg/K.

  Returns:
    double: density following Eq. (4.2.8) in Emanuel (1994).
%}
  global RD
  rho = p / RD ./ eval_Tv(T, q, ql);