function theta_e = eval_theta_e(T, q, ql, qi, p, P0)
%{
  Calculate equivalent potential temperature.

  Units: K.

  Following Eq. (2.67) in Stevens & Siebesma (2020)
   which, with ice correction, is equivalent to Eq. (4.5.11) in Emanuel (1994).

  Args:
    T (double): air temperature in K.
    q (double): specific humidity of water vapor in kg/kg.
    ql (double): specific liquid water content in kg/kg.
    qi (double): specific liquid water content in kg/kg.
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
    global RD RV CPD CL LF EPS
    if nargin==5
      P0 = 1e5;
    end
    T0 = 273.1636783445389; % Triple point = freezing point; for water-ice saturation consistency
    LV = lv(T);
    
    qt = q+ql+qi;
    e = p.*q./( q+EPS*(1-qt) );
    es = eval_es(T).*(T>=T0)+eval_esi(T).*(T<T0);
    
    RE = (1-qt)*RD;
    R = RE + q*RV;
    
    CPL = CPD + qt*(CL-CPD);
    
    chi = RE./CPL;
    gamma = RV*q ./ CPL;
    gammi = (CL-ci(T0)) .* qi ./CPL;
    
    theta_e = T .* (P0*R./p./RE).^chi ...
                .* (es./e).^gamma ...
                .* (T0./T).^gammi ...
                .* exp( (q.*LV./T - LF*qi/T0)./CPL );
              
%% Below: Eq. (4.5.11) from Emanuel (1994); equivalent to Eq. (2.67) in Stevens & Siebesma (2020)  
%     r = q ./ (1-qt);
%     ri = qi ./ (1-qt);
%     rt = qt ./ (1-qt);
%     
%     e = p.*q./( q+EPS*(1-qt) );
%     es = eval_es(T).*(T>=T0)+eval_esi(T).*(T<T0);
%     rh = e./es;
%     pd = p-e;
%     
%     CP = CPD + rt*CL;
%     chi = RD./CP;
%     gamma = -r.*RV./CP;
%     
%     theta_e = T .* (P0./pd).^chi ...
%                 .* rh.^gamma ...
%                 .* (T/273.16).^( -ri.*(CL-ci(T))./CP ) ...
%                 .* exp( (LV.*r./T - LF*ri/273.16)./CP );