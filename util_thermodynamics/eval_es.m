function es = eval_es(T)
%{
  Calculate saturation vapor pressure with respect to liquid water.

  Units: Pa.

  Args:
    T (double): air temperature in K.

  Returns:
    double: saturation vapor pressure of water (liquid) using Eq. (4.4.13) in Emanual (1994) 
            Error within 0.006% in 0<=T<=40-deg-C; 
                           0.3% for -30-deg-C<=T;
                           0.7% for -40-deg-C<=T.
%}
es = 1e2 * exp( 53.67957 - 6743.769./T - 4.8451*log(T) );
% d(es)/dT = (6743.769./T-4.8451).*es./T