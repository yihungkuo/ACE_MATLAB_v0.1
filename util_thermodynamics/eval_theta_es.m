function theta_es = eval_theta_es(T, p, P0)
%{
  Calculate saturation equivalent potential temperature.

  Units: K.

  Following Eq. (2.67) in Stevens & Siebesma (2020)
   which, with ice correction, is equivalent to Eq. (4.5.11) in Emanuel (1994).

  Args:
    T (double): air temperature in K.
    p (double): air pressure in Pa.
    RD (constant): gas constant of dry air in J/Kg/K.
    RV (constant): gas constant of water vapor in J/Kg/K.
    CPD (constant): heat capacity at constant pressure for dry air in J/Kg/K.
    CL (constant): heat capacity of liquid water (above freezing) in J/Kg/K.
    EPS (constant): ratio of gas constant of dry air to gas constant of water vapor.
    P0 (constant): reference pressure (1,000 hPa) in Pa.

  Returns:
    double: equivalent potential temperature.
	%}
    global RD RV CPD CL EPS
    if nargin==2
      P0 = 1e5;
    end
    T0 = 273.1636783445389; % Triple point = freezing point; for water-ice saturation consistency
    LV = lv(T);
    
    q = eval_qsat(T,p);
    qt = q;
    e = p.*q./( q+EPS*(1-qt) );
    es = eval_es(T).*(T>=T0)+eval_esi(T).*(T<T0);
    
    RE = (1-qt)*RD;
    R = RE + q*RV;
    
    CPL = CPD + qt*(CL-CPD);
    
    chi = RE./CPL;
    gamma = RV*q ./ CPL;
    
    theta_es = T .* (P0*R./p./RE).^chi ...
                .* (es./e).^gamma ...
                .* exp( (q.*LV./T)./CPL );  