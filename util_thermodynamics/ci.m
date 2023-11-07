function ciT = ci(T)
%{
	Calculate heat capacity of ice.

  Units: J/Kg/K.
    
  Args:
    T (double): air temperature in K.
        
	Returns:
    double: heat capacity of ice at temperature T following Appendix 2 (p. 566-567) of Emanuel (1994).
%}
	ciT = 2106 + 7.3*(T-273.15);