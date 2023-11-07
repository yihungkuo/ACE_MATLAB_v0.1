function [T, q, ql, qi] = inv_T_from_theta_il(theta_il, qt, p, ABS_TOL, P0)
% function [T, q, ql, qi, step] = inv_T_from_theta_il(theta_il, qt, p, ABS_TOL, P0)
%{
  Invert air temperature and specific humidities, 
  given liquid-ice water potential temperature, total specific water
  content, pressure.

  Units: K.

  Args:
    theta_il (double): liquid-ice water potential temperature in K.
    qt (double): specific total water content in kg/kg.
    p (double): air pressure in Pa.
    RD (constant): gas constant of dry air in J/Kg/K.
    RV (constant): gas constant of water vapor in J/Kg/K.
    CPD (constant): heat capacity at constant pressure for dry air in J/Kg/K.
    CPV (constant): heat capacity at constant volume for dry air in J/Kg/K.
    CPVMCL (constant): heat capacity of liquid water minus CPV.
    LS (constant): latent heat of sublimation in J/Kg/K (-100<=T<= 0-deg-C)
    EPS (constant): ratio of gas constant of dry air to gas constant of water vapor.
    P0 (constant): reference pressure (1,000 hPa) in Pa.
    ABS_TOL (constant): absolute tolerance (units: K) for stopping iteration.

  Returns:
    double: air temperature.
    double: specific humidities.
    step: for debugging
%}
%% Goal: solve (T,q,ql,qi) given (theta_il,p,qt)
%    F(T,...) = c1 * T * R^chi / qv^e1 * exp( E(T,ql,qi) ) - theta_il == 0  ...Eq. (1)
global CPD CPV RD RV CPVMCL LS EPS

if nargin==3
  ABS_TOL = 1e-3; %<--Apply only for Newton iteration: NOT for ice-liquid coexisting
  P0 = 1e5;
elseif nargin==4
  P0 = 1e5;
end

T0 = 273.1636783445389; % Triple-point temperature; for water-ice saturation consistency

CPL = CPD + (CPV-CPD)*qt;
RL = RD*(1-qt) + RV*qt;
chi = RL./CPL;

%% Step 1: solve assuming no condensate & check consistency
% If qc>0, T would be lower than true value by a factor~exp(E(T,ql,qi))
PI = (p./P0).^chi; % RL=R for condensate-free
T = theta_il.*PI;
q = qt;
ql = 0;
qi = 0;
% step = 1;

% Check no-condensate consistency: if condensate, qt>qs
if T>=T0
  qsat = eval_qs(T,p,qt);
else
  qsat = eval_qsi(T,p,qt);
end

if q>qsat % If condensate ==> Step 2
  %% Step 2: solve assuming liquid/no-ice & check consistency
  %  From Eq. (1): G(T) = T * R^chi / qs(T)^e1 * exp( -E(T,qt-qs(T),qi=0) ) - theta_il/c1 = 0  ...Eq. (2)
  gamma = RV*qt./CPL;
  G0 = theta_il .* PI .* RL.^chi ./ qt.^gamma; % G0=theta_il/c1
  
  % Newton would fail if initial T too small (since dG/dT depends ~exponential on T)
  % Set initial T s.t. qt=qs(T)~EPS*es(T)/p with Tetens: es(T)~610.78*exp(17.27*(T-273.16)/(T-35.86))
  LOG = log(p.*qt) - log(610.78*EPS);
  T = 237.3*LOG./(17.27-LOG) + 273.16;
  
  dT = 9999;
  while abs(dT)>ABS_TOL
    qs = eval_qs(T,p,qt);
    dTqs = eval_dTqs(T,p,qt);
    ql = qt-qs;
    LV = lv(T);
    R = RD*(1-qt) + RV*qs;
    % The following 3 lines for dG/dT verified!
    EXP = exp(-ql.*LV./CPL./T);
    dGdT = EXP .* R.^chi ./ CPL ./ qs.^gamma ...
               .*( dTqs.*( LV + RV*(RL./R - qt./qs).*T ) + CPL + ql.*(LV./T+CPVMCL) );
    G = T./qs.^gamma .* EXP .* R.^chi - G0;
    % Newton: T(n+1) = T(n) - G(Tn)/G'(Tn)
    dT = -G./dGdT;
    T = T + dT;
  end%convergence
  q = eval_qs(T,p,qt);
  ql = qt-q;
  qi = 0;
  if ql<0
    q = qt;
    ql = 0;
  end
%   step = 2;

  % Check liquid/no-ice consistency: if ice exists (so T<=T0) --> T<T0
  if T<T0 % If ice ==> Step 3
    %% Step 3: solve assuming ice/no-liquid & check consistency
    % Tetens for ice: esi(T)~610.78*exp(21.875*(T-273.16)/(T-7.66))
    LOG = log(p.*qt) - log(610.78*EPS);
    T = 265.5*LOG./(21.875-LOG) + 273.16;
    
    dT = 9999;
    while abs(dT)>ABS_TOL
      qsi = eval_qsi(T,p,qt);
      dTqsi = eval_dTqsi(T,p,qt);
      qi = qt-qsi;
      R = RD*(1-qt) + RV*qsi;
      % The following 3 lines for dG/dT verified!
      EXP = exp(-LS*qi./CPL./T);
      dGdT = EXP .* R.^chi ./ CPL ./ qsi.^gamma ...
                 .*( dTqsi.*( LS + RV*(RL./R - qt./qsi).*T ) + CPL + LS*qi./T );
      G = T./qsi.^gamma .* EXP .* R.^chi - G0;
      % Newton: T(n+1) = T(n) - G(Tn)/G'(Tn)
      dT = -G./dGdT;
      T = T + dT;
    end%convergence
    q = eval_qsi(T,p,qt);
    qi = qt-q;
    ql = 0;
    if qi<0
      q = qt;
      qi = 0;
    end
%     step = 3;

    % Check ice/no-liquid consistency: if liquid exists (so T>=T0), T>T0 (and qc<0 for some cases)
    if T>T0 % If ice-liquid coexisting ==> Step 4
      %% Step 4: solve assuming T=T0, ice-liquid coexisting
      T = T0;
      q = (eval_qs(T0,p,qt) + eval_qsi(T0,p,qt))/2; %<--qs(T0)=qsi(T0) in theory but... ==> use mid-point value
      R = RD*(1-qt) + RV*q;
      % Just find out the partition (ql,qi) s.t. ql+qi=qt-qs
      qc = qt-q;
      ql = ( LS*qc - CPL.*T.*log((T0./G0) ./ q.^gamma .* R.^chi) )./(LS-lv(T0));
      qi = qc-ql;
      if qc<0 %if discrepancy arises from qs(T0)>qsi(T0) [because of formulae for es(T) & esi(T)]
        q = qt;
        ql = 0;
        qi = 0;   
      elseif ql<0 %if liquid-free
        ql = 0;
        qi = qc;
      elseif qi<0 %if ice-free
        qi = 0;
        ql = qc;
      end
%       step = 4;

    end%if ice/no-liquid
  end%if liquid/no-ice
end%if condensate