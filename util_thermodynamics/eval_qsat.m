function qsat = eval_qsat(T, p, qt, T0)
%{
  Calculate saturation specific humidity.

  Depending on temperature, this subroutine 
  calls eval_qs (T>273.15K) & eval_qsi (T<=273.15K).

  Units: kg/kg (unitless).

  Args:
    T (double): air temperature in K.
    p (double): air pressure in Pa.
    qt (double): total water mixing ratio in kg/kg.
    T0 (double): freezing point in K.

  Returns:
    double: saturation specific humidity.
%}
  if nargin==2  %assuming condensate-free & freezing at triple point
    T0 = 273.1636783445389; %for consistent water and ice saturation at triple point
    qsat = eval_qs(T, p);
    qsat(T<T0) = eval_qsi(T(T<T0), p(T<T0));
  elseif nargin==3
    T0 = 273.16;
    qsat = eval_qs(T, p, qt);
    qsat(T<T0) = eval_qsi(T(T<T0), p(T<T0), qt(T<T0));
  else
    qsat = eval_qs(T, p, qt);
    qsat(T<T0) = eval_qsi(T(T<T0), p(T<T0), qt(T<T0));
  end
  
  