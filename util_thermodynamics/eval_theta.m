function theta = eval_theta(T, p, P0)
%{
  Calculate (dry) potential temperature.

  Units: K.

  Args:
    T (double): air temperature in K.
    p (double): air pressure in Pa.
    RD (constant): gas constant of dry air in J/Kg/K.
    CPD (constant): heat capacity at constant pressure for dry air in J/Kg/K.
    P0 (constant): reference pressure (1,000 hPa) in Pa.

  Returns:
    double: potential temperature.
	%}
    global RD CPD
    if nargin==2
      P0 = 1e5;
    end
    
    theta = T .* (P0./p).^(RD/CPD);