function esi = eval_esi(T)
%{
  Calculate ice saturation vapor pressure.

  Units: Pa.

  Accuracy within 0.14% for -80<=T<=0-deg-C.

  Args:
    T (double): air temperature in K.

  Returns:
    double: ice saturation vapor pressure following Eq. (4.4.15) in Emanuel (1994).
            Error within 0.14% for -80<=T<=0-deg-C.
%}   
	esi = 100 * exp( 23.33086 - 6111.72784./T + 0.15215*log(T) );