function lvT = lv(T)
%{
  Calculate latent heat of vaporization

  Units: J/Kg.

  Args:
    T (double): air temperature in K.
  
  Returns:
    double: latent heat of vaporization at temperature T following Appendix 2 (p. 566-567) of Emanuel (1994).
    
%}
  global LV0 CPVMCL
  lvT = LV0 - CPVMCL*(T-273.15);