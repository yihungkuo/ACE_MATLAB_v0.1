% =========================================================================
% ace_v01.m
%
%   Anelastic Convective Entity (ACE) Model - MATLAB Implementation Version 0.1
%   06-Nov-2023
%   Model development: Yi-Hung Kuo & J. David Neelin
%   Code development: Yi-Hung Kuo
%   yhkuo@princeton.edu (yhkuo@atmos.ucla.edu)
%   CIMES, Princeton University
%
%   This is the driver script for the ACE_MATLAB package prepared for and
%   released with the submission of Kuo & Neelin (2023) to JAS. 
%
%   Utility subroutines used by this driver script can be found under
%    ACE_MATLAB_v1.0/util_weno/
%    ACE_MATLAB_v1.0/util_thermodynamics/
%   ARMBE soundings under
%    ACE_MATLAB_v1.0/data/soundings/  <-----NOT PART OF THE ACE_MATLAB PACKAGE!
%   Pre-computed nonlocal bases under
%    ACE_MATLAB_v1.0/data/basis/
%
%   Reference: 
%    Kuo, Y.-H., and J. D. Neelin, 2022: Conditions for convective deep inflow. 
%      Geophys. Res. Lett., 49, e2022GL100 552, doi:10.1029/2022GL100552.
%    Kuo, Y.-H., and J. D. Neelin, 2023: Anelastic Convective Entities:
%      formulation and properties. J. Atmos. Sci., submitted. Preprint
%      available at https://drive.google.com/file/d/1e0_n11Eop3eqQGu_-Vtl7PI2V7KjTI33/view?usp=sharing.
% =========================================================================

%% Directory/file paths
dir_ace = pwd;

% Subroutines
dir_util_weno = [dir_ace '/util_weno/'];
% dir_util_plume = [dir_ace '/util_plume/'];
dir_util_thermodynamics = [dir_ace '/util_thermodynamics/'];
addpath(dir_util_weno,dir_util_thermodynamics)

% Pre-computed anelastic basis for massflux
dir_basis = [dir_ace '/data/basis/'];

% ARM best-estimate soundings (3 files for 3 sites)
dir_armbe = [dir_ace '/data/soundings/'];
armname = {'Manus';'Nauru';'GoAmazon'};
armfnlist = {'manus_3hr_ARMBEsonde_profiles.nc';
             'nauru_3hr_ARMBEsonde_profiles.nc';
             'maoarmbeatmM1.c1.20140101.003000.nc'};

% Output folder
ace_output_dir = [dir_ace '/ace_output/'];

%% Declare global for in-script functions
global theta_e_mf qt_mf Tv_mf p_mf rho_mf rho_mfh MB Z zm z p g k dz Nzm 
global eps_tur eps_sl eddiffu eddiwin qc_ramp dqc tau_pr EBF EBT

%% Set up parameters
plotting = 1; %default colormap=turbo (for R202x or later)

% ARM sounding option
arm = 3; %1=Manus/2=Nauru/3=GOAmazon
armfn = armfnlist{arm};
pidx = 314; % Profile index for environment

% Output filename
ace_output_fn = [armname{arm} '_pidx=' num2str(pidx) '.mat']; % Feel free to append

% Horizontal entity diameter
D = 5; % km 
acylin = 0;
if acylin==0 % Loading horizontally cylindrical basis 
  pcbsfn = ['PRECOMPUTED_BASIS_D=' num2str(D) 'KM.mat'];
else % Acylindrical (elliptical) basis 
  pcbsfn = ['PRECOMPUTED_BASIS_ACYLIN_D~' num2str(D) 'KM.mat'];
end

% Verital grids setup
dz = 50; % m; 5hPa~50m(surface)/250m(aloft)
ztop = 40e3; % m
z = [dz/2:dz:ztop-dz/2]'; % z-grid for plume computation
zm = [0:dz:ztop]'; % z-grid for mass-flux for entrainment rate
p = [1000:-5:100]'*1e2; % 5-hPa grid

% For precip sinks
qc_ramp = 0.5e-3; % kg/kg; ramping for precip when qc>qc_ramp
dqc = 20; % kg/kg; transition for a smooth onset of precip ramp function
tau_pr = 5*60; % sec; relaxation timescale for precip

% For ODE solver
conti = 0; % Restart or continue a run
% Initiate for continuing runs
if conti==1
  Xminus1 = X(:,1:end-1);
  tminus1 = tspan(1:end-1);
  X0 = X(:,end);
  t0 = tspan(end); % t0=starting time; T0=triple point temperature.
  disp('Continue Previous Run...')
else
  disp('Start Anew!')
end
% DO NOT MODIFY THE IF STATEMENTS
tspan = [0:0.25:90]*60; % sec
Nt = length(tspan);
% Update tspan for continuing runs
if conti==1
  tspan = tspan + t0;
end
% DO NOT MODIFY THE IF STATEMENTS
MaxStep = 1; % MaxStep for ode solver

% Initial mass flux = MF0*sin(pi*(z-zb)/(zt-zb)) for zb<z<zt (zeros elsewhere)
MF0 = 5; % kg/m^2/s
zt = 2e3; % km; top of initial MF
zb = 0e3; % km; Base of inital MF

% Initial thermal bubble in zm(BL); BL=[] for no bubble
BL = []; % BL=bubble levels
dTbl = 0; % T=T+dTbl in BL
dqsat = 1; %1=saturate zm(BL); 0=doing nothing
dqcbl = 0;%10e-3; % kg/kg; qc=dqcbl; ignored if dqsat~=1
dqbl = 0.8; % q=dqbl*q in BL; ignored if dqsat==1

% External buoyancy forcing (mimicking large-scale forcing)
EB = 0.02e-2; % m/s^2; magnitude of external buoyancy
zebb = 0; % Base of external buoyancy  layer
zebt = 2e3; % Top of external buoyancy  layer
EBF = EB * normpdf(zm,(zebt-zebb)/2,(zebt-zebb)/7) / normpdf(0,0,(zebt-zebb)/7);
EBT = 0*60; % sec; EBF turned on for t<EBT

% Turbulent entrainment rate
nju = 0.09; % Morrison (2017)
eps_tur = 4*nju/D/1e3; % Note: D in km

% Eddy diffusivity substituting additional dynamic-induced pressure perturbation terms
eddiffu = D*1e3/4/pi/dz;
eddiwin = ceil(D*1e3/2/dz); % Sliding window half-width for movmax

% Spounge layer with Newtonian damping at domain top ztop
SLD = 1e3;%3e3; % e-folding depth in m
SLTau = 60; % SL relaxation timescale in s
eps_sl = exp( (zm-zm(end))/SLD )/SLTau;

% WENO scheme (reconstruction uses 2k-1 points)
k = 3;

% Other pre-defined constants
IMPORT_CONSTANTS
Nzm = length(zm);

%% Loading pre-computed massflux basis
load([dir_basis pcbsfn])

%% Pre-processing sounding data
[plev, pr, ta, hus, date] = LOAD_SOUNDINGS([dir_armbe armfn], arm); %<---GoAmazon midnight hardwired by hh==5

% Interpolate sounding values to 5-hPa p-grid (_p)
[T_p, q_p, rho_p, z_p] = INTERP_P_GRID(plev, ta, hus);

% Interpolate sounding values to 50-m z-grid (_env) & zm-grid (_zm)
[T_env, q_env, Tv_env, theta_e_env, rho_env, p_z,...
 T_zm, q_zm, Tv_zm, theta_e_zm, rho_zm, p_zm,...
 theta_zm, theta_es_zm, rh_zm] = INTERP_Z_GRID(T_p, q_p, z_p);

%% Setup for time-integration
day = date(pidx); % Can run multiple pidx cases (by adding a for-loop here)

% _mf(h) = properties of mass flux; will be used as constant reference during time-integration
theta_e_mf = theta_e_zm(:,pidx); 
qt_mf = q_zm(:,pidx);
Tv_mf = Tv_zm(:,pidx);
p_mf = p_zm(:,pidx);
rho_mf = rho_zm(:,pidx); 
rho_mfh = interp1(zm,rho_mf,zm+dz/2,'pchip');

% Mass-flux basis MB(Z,Bi)=PMB*rho_0(Z)
rho_Z = rho_mf(1:4:end);
MB = PMB.*rho_Z;

% Initiate MF(t=0)
MF = MF0*sin( pi*(zm-zb)/(zt-zb) );
MF(zm>zt | zm<zb) = 0;

% Initiate for ODE solver
if conti == 0
  X0 = [theta_e_mf,qt_mf,MF,zeros(Nzm,2)]; %<---variables defined on zm-grid
  t0 = 0;
  % Adding BL/thermal bubble
  if ~isempty(BL)
    if dqsat==1
      X0(BL,2) = eval_qsat(T_zm(BL,pidx)+dTbl,p_mf(BL)) .* (1-dqcbl);
    else
      X0(BL,2) = q_zm(BL,pidx)*dqbl;
    end
    X0(BL,1)=eval_theta_e(T_zm(BL,pidx)+dTbl,X0(BL,2),(dqsat==1)*dqcbl*((T_zm(BL,pidx)+dTbl)>=T0),(dqsat==1)*dqcbl*((T_zm(BL,pidx)+dTbl)<T0),p_mf(BL));
  end
end

%% Time integration
opt = odeset('NonNegative',Nzm+(1:Nzm),'MaxStep',MaxStep);

%  X = [theta_e_mf, qt_mf, MF, dtheta_e_pr_accum, pr_accum];
[t,X] = ode45(@odeFun,tspan,X0(:),opt);
X = X';
if conti==1
  X = [Xminus1 X];
  tspan = [tminus1 tspan];
  Nt = length(tspan);
end

% Derived variables evaluated via post-processing
[theta_e, qt, mf, theta, RH, B,...
 qsat, qc, ql, qi, Tv, T, dtheta_e_pr_accum, pr_accum] = POST_PROCESS(X);

% Save output
mkdir(ace_output_dir)
save([ace_output_dir ace_output_fn])

%% Plotting & Saving
if plotting

  % Equivalent potential temperature
  figure('Units','normalized','Position',[0 0 1 1])
  contourf(tspan/60,zm/1e3,theta_e,[330:0.5:360],'linestyle','none')
  hold all
  colormap(turbo)
  contour(tspan/60,zm/1e3,qc*1e3,[0.01 0.01],'color','w','linewidth',2)
  ylim([0,14])
  grid on
  clb = colorbar;
  clb.Ticks = [335:3:350];
  clim([335 350])
  set(gca,'fontsize',30)
  xlabel('Time (min)')
  ylabel('Height (km)')
  yticks(0:1:20)
  yticklabels({'0';'';'2';'';'4';'';'6';'';'8';'';'10';'';'12';'';'14'})
  ax = gca;
  ax.GridAlpha = 0.5;
  title(['ACE {\it \theta_e (K)} - ' armname{arm} ' ' datestr(day) ' (UTC)'])

  % Specific mass of condensate
  figure('Units','normalized','Position',[0 0 1 1])
  contourf(tspan/60,zm/1e3,qc*1e3,[0.01,0.05:0.05:25],'linestyle','none')
  hold all
  colormap(turbo)
  hold all
  ylim([0,14])
  grid on
  clb = colorbar;
  clb.Ticks = [0.01,1:1:5];
  clim([0.01 5])
  set(gca,'fontsize',25)
  xlabel('Time (min)')
  ylabel('Height (km)')
  yticks(0:1:20)
  yticklabels({'0';'';'2';'';'4';'';'6';'';'8';'';'10';'';'12';'';'14'})
  ax = gca;
  ax.GridAlpha = 0.5;
  title(['ACE {\it q_c (g/kg)} - ' armname{arm} ' ' datestr(day) ' (UTC)'])

  % Mass flux
  figure('Units','normalized','Position',[0 0 1 1])
  contourf(tspan/60,zm/1e3,mf,[-4:0.1:12],'linestyle','none')
  hold all
  colormap(turbo)
  contour(tspan/60,zm/1e3,qc*1e3,[0.01 0.01],'color','w','linewidth',2)
  ylim([0,14])
  grid on
  clb = colorbar;
  clb.Ticks = [-3:1:8];
  clim([-3 8])
  set(gca,'fontsize',20)
  xlabel('Time (min)')
  ylabel('Height (km)')
  yticks(0:1:20)
  yticklabels({'0';'';'2';'';'4';'';'6';'';'8';'';'10';'';'12';'';'14'})
  ax = gca;
  ax.GridAlpha = 0.5;
  title(['ACE {\it \rho_0w (kg/m^2/s)} - ' armname{arm} ' ' datestr(day) ' (UTC)'])

  % d(Mass flux)/dz = dynamic entrainment/detrainment
  figure('Units','normalized','Position',[0 0 1 1])
  contourf(tspan/60,zm/1e3,dMFdz,[-5:0.1:5]*1e-3,'linestyle','none')
  colormap(turbo)
  hold all
  contour(tspan/60,zm/1e3,qc*1e3,[0.01 0.01],'color','w','linewidth',2)
  ylim([0,14])
  grid on
  clb = colorbar;
  clb.Ticks = [-3.5:0.5:2.5]*1e-3;
  clim([-3.5 2.5]*1e-3)
  set(gca,'fontsize',20)
  xlabel('Time (min)')
  ylabel('Height (km)')
  yticks(0:1:20)
  yticklabels({'0';'';'2';'';'4';'';'6';'';'8';'';'10';'';'12';'';'14'})
  ax = gca;
  ax.GridAlpha = 0.5;
  title(['ACE {\it \partial_z(\rho_0w) (kg/m^3/s)} - ' armname{arm} ' ' datestr(day) ' (UTC)'])

  % Precipitation
  figure
  colormap(jet)
  subplot(1,3,1)
  contourf(tspan/60,zm/1e3,dtheta_e_pr_accum,'LineStyle','none')
  grid on
  title('\Delta\theta_e (K)')
  xlabel('Time (min)')
  ylabel('Height (km)')
  colorbar
  set(gca,'fontsize',20)
  subplot(1,3,2)
  contourf(tspan/60,zm/1e3,pr_accum,'LineStyle','none')
  grid on
  title('Accum. Precip. (mm)')
  xlabel('Time (min)')
  ylabel('Height (km)')
  colorbar
  set(gca,'fontsize',20)
  subplot(1,3,3)
  plot(tspan/60,sum(pr_accum)')
  grid on
  xlabel('Time (min)')
  ylabel('Accum. Precip. (mm)')

end%plotting

%end %<--- When running multiple pidx cases, for-loop ends here

%% Subroutines
%%   implemented
%%     below:
function dXdt = odeFun(t,X)
global theta_e_mf qt_mf rho_mf rho_mfh k dz Nzm eps_tur eps_sl eddiffu eddiwin EBF EBT
  disp(t)
  X = reshape(X,[Nzm,5]); % X=[theta_e,qt,MF,dtheta_e_pr_accum,pr_accum]
  Xext = [ [ones(k-1,1)*X(1,1);X(:,1)],...
           [ones(k-1,1)*X(1,2);X(:,2)],...
           [-X(k:-1:2,3);X(:,3)] ];
  Nzext = Nzm+k-1; %=length(Xext)/3
  dXextdz = zeros(Nzext,3);

  % Compute MFh(zm+dz/2) from MF(zm) as the mean of up- and down-wind approximations
  LMFh = zeros(Nzm,1); % MFh from more left stencil points
  RMFh = zeros(Nzm,1); % MFh from more right stencil points
  for idz=k+1:Nzext-k+1
    LMFh(idz-k) = reconstruction_weno(k, Xext([idz-k:idz+k-2],3) );
    RMFh(idz-k) = reconstruction_weno(k, Xext(flip([idz-k+1:idz+k-1]),3) );
  end
  MFh = (LMFh + RMFh)/2;
  % w at zm+dz/2 with up-/down-wind following Godunov Flux
  wh = ( LMFh.*(MFh>=0) + RMFh.*(MFh<0) )./rho_mfh;
  for idz=k+1:Nzext-k
    % d(theta_e)/dz
    dXextdz(idz,1)...
        = ( reconstruction_weno(k, Xext( [idz-k+1:idz+k-1]*(MFh(idz-k+1)>=0)+flip([idz-k+2:idz+k])*(MFh(idz-k+1)<0) ,1) ) ...
           - reconstruction_weno(k, Xext( [idz-k:idz+k-2]*(MFh(idz-k)>=0)+flip([idz-k+1:idz+k-1])*(MFh(idz-k)<0) ,1) ) )/dz;
    % d(qt)/dz
    dXextdz(idz,2)...
        = ( reconstruction_weno(k, Xext( ([idz-k+1:idz+k-1])*(MFh(idz-k+1)>=0)+flip([idz-k+2:idz+k])*(MFh(idz-k+1)<0) ,2) ) ...
           - reconstruction_weno(k, Xext( ([idz-k:idz+k-2])*(MFh(idz-k)>=0)+flip([idz-k+1:idz+k-1])*(MFh(idz-k)<0) ,2) ) )/dz;
    % d(MF=rho_0*w)/dz
    dXextdz(idz,3)...
        = ( reconstruction_weno(k, Xext( ([idz-k+1:idz+k-1])*(MFh(idz-k+1)>=0)+flip([idz-k+2:idz+k])*(MFh(idz-k+1)<0) ,3) ) ...
           - reconstruction_weno(k, Xext( ([idz-k:idz+k-2])*(MFh(idz-k)>=0)+flip([idz-k+1:idz+k-1])*(MFh(idz-k)<0) ,3) ) )/dz;
  end
  % First zero at zm=0 (dXdz~=0 but dXdt=0 at the end of computation)
  % zeros(k,1) for domain top<---to be improved!!!!!
  dXdz = [ [0;dXextdz(k+1:Nzext-k,1);zeros(k,1)],...
           [0;dXextdz(k+1:Nzext-k,2);zeros(k,1)],...
           [0;dXextdz(k+1:Nzext-k,3);zeros(k,1)] ];
  % E: dynamic & turbulent entrainment contributions
  eps_dyn = (dXdz(:,3)>0).*dXdz(:,3);
  E1 = ( theta_e_mf-X(:,1) ) .*( eps_tur.*rho_mf + eps_sl + eps_dyn );
  E2 = ( qt_mf-X(:,2) ).*( eps_tur.*rho_mf + eps_sl + eps_dyn );
  % E3 includes advection of MF (as an input for nonLocalResponse)
  E3 = -(wh(2:end).^2-wh(1:end-1).^2)/(2*dz)...
        -X(2:end,3)./(rho_mf(2:end).^2).*( eps_tur.*rho_mf(2:end) + eps_sl(2:end) + eps_dyn(2:end) );
  % Eddy Diffusion for Momentum
  ED = abs(dXdz(:,3))*( eddiffu*(max(X(:,3)./rho_mf)-min(X(:,3)./rho_mf))/max(abs(dXdz(:,3))+eps)/dz );
  ED = movmax(ED,[1 1]*eddiwin);
  ED = conv(ED,ones(2*eddiwin+1,1)/(2*eddiwin+1),'same');
  ED = ED .* ( Xext(k-1:end-1,3)-2*X(:,3)+[Xext(k+1:end,3);-Xext(end,3)] );
  % Ivert (T,q,ql,qi) for buoyancy and precipitation
  [T, Q, QL, QI] = INV_THETA_E(X(:,1), X(:,2));
  % Precipitation sinks
  [dtheta_e, dqt] = PRECIP_SINKS(T, Q, QL, QI);
  dXdt = [ ( E1 - dXdz(:,1).*X(:,3) )./rho_mf + dtheta_e,...
           ( E2 - dXdz(:,2).*X(:,3) )./rho_mf + dqt,...
           NONLOCAL_RESPONSE(T, Q, QL+QI, [0;E3]) + ED + EBF*(t<EBT),...
           dtheta_e, -dqt];
  dXdt([1,end],:) = 0; % At zm=0 & domain top
  dXdt = dXdt(:);
end

function dMFdt = NONLOCAL_RESPONSE(T, q, qc, E3)
global Tv_mf MB Z zm z g
  Tv = eval_Tv(T, q, qc);
  B = g.*(Tv-Tv_mf)./Tv_mf;
  % dMFdt = L(B;D)<---averaged MF tendency in r<D/2
  dMFdt = interp1(Z,MB*mean(reshape(interp1(zm,B+E3,z),4,[]))',zm,'pchip');
  dMFdt(1) = 0;
end

function [dtheta_e, dqt] = PRECIP_SINKS(T, q, ql, qi)
global p_mf qc_ramp dqc tau_pr
  qc = ql+qi;
  dqt = -(qc_ramp/dqc/tau_pr) * log( 1+exp( (dqc/qc_ramp) * (qc-qc_ramp) ) );
  % Finite-difference estimate for d(theta_e)/dqc
  dql = 1e-8*ql./qc;
  dqi = 1e-8*qi./qc;
  dtheta_e = dqt.* ( eval_theta_e(T,q,ql+dql,qi+dqi,p_mf) - eval_theta_e(T,q,ql,qi,p_mf) )/1e-8;
  % Check qc>0
  dqt(qc<=0) = 0;
  dtheta_e(qc<=0) = 0;
end

function [T, q, ql, qi] = INV_THETA_E(theta_e, qt)
global p_mf Nzm
  T = zeros(Nzm,1);
  q = zeros(Nzm,1);
  ql = zeros(Nzm,1);
  qi = zeros(Nzm,1);
  for idz=1:Nzm
    [T(idz), q(idz), ql(idz), qi(idz)] = inv_T_from_theta_e(theta_e(idz), qt(idz), p_mf(idz));
  end
  % Double check values after inversion
  T = real(T);
  q = real(q);
  ql = real(ql);
  qi = real(qi);
end

function [theta_e, qt, mf, theta, RH, B,... 
          qsat, qc, ql, qi, Tv, T, dtheta_e_pr_accum, pr_accum] = POST_PROCESS(X)
global Tv_mf p_mf Nzm g
  theta_e = X(1:Nzm,:);
  qt = X(Nzm+1:2*Nzm,:);
  mf = X(2*Nzm+1:3*Nzm,:);
  dtheta_e_pr_accum = X(3*Nzm+1:4*Nzm,:);
  pr_accum = X(4*Nzm+1:end,:);
  Nt = size(X,2);
  theta = zeros(Nzm,Nt);
  qsat = zeros(Nzm,Nt);
  ql = zeros(Nzm,Nt);
  qi = zeros(Nzm,Nt);
  T = zeros(Nzm,Nt);
  for idt=1:Nt
    for idz=1:Nzm
      [T(idz,idt), q(idz,idt), ql(idz,idt), qi(idz,idt)]...
          = inv_T_from_theta_e(theta_e(idz,idt), qt(idz,idt), p_mf(idz));
    end
    theta(:,idt) = eval_theta(T(:,idt),p_mf);
    qsat(:,idt) = eval_qsat(T(:,idt),p_mf,qt(:,idt));
  end
  Tv = eval_Tv(T,q,ql+qi);
  RH = q./qsat;
  B = g.*(Tv-Tv_mf)./Tv_mf;
  qc = (qt-qsat).*(qt>qsat);
end

function [T_env, q_env, Tv_env, theta_e_env, rho_env, p_z,...
          T_zm, q_zm, Tv_zm, theta_e_zm, rho_zm, p_zm, theta_zm,...
          theta_es_zm, rh_zm] = INTERP_Z_GRID(T_p, q_p, z_p)
global RD g dz z zm p
  % Interpolate sounding values to 50-m z-grid (_env) & zm-grid (_zm)
  T_env = nan([length(z),size(T_p,2)]);
  q_env = nan([length(z),size(T_p,2)]);
  Tv_env = nan([length(z),size(T_p,2)]);
  p_z = nan([length(z),size(T_p,2)]);
  T_zm = nan([length(zm),size(T_p,2)]);
  q_zm = nan([length(zm),size(T_p,2)]);
  Tv_zm = nan([length(zm),size(T_p,2)]);
  p_zm = nan([length(zm),size(T_p,2)]);
  for pidx=1:size(T_p,2)
    T_env(:,pidx) = interp1(z_p(:,pidx),T_p(:,pidx),z,'pchip',T_p(end,pidx));
    q_env(:,pidx) = interp1(z_p(:,pidx),q_p(:,pidx),z,'pchip',q_p(end,pidx));
    q_env(q_env(:,pidx)<0,pidx) = 0; % Fix negative q_env due to interp1
    Tv_env(:,pidx) = eval_Tv(T_env(:,pidx),q_env(:,pidx));
    p_z(:,pidx) = interp1(z_p(:,pidx),p,z,'pchip');
    midz = sum(z<=max(z_p(:,pidx)));
    for idz=midz+1:length(z)
      p_z(idz,pidx) = p_z(idz-1,pidx)*(1-g*dz/2/RD/Tv_env(idz-1,pidx))/(1+g*dz/2/RD/Tv_env(idz-1,pidx));
    end
    T_zm(:,pidx) = interp1(z_p(:,pidx),T_p(:,pidx),zm,'pchip',T_p(end,pidx));
    q_zm(:,pidx) = interp1(z_p(:,pidx),q_p(:,pidx),zm,'pchip',q_p(end,pidx));
    q_zm(q_zm(:,pidx)<0,pidx) = 0; % Fix negative q_zm due to interp1
    Tv_zm(:,pidx) = eval_Tv(T_zm(:,pidx),q_zm(:,pidx));
    p_zm(:,pidx) = interp1(z_p(:,pidx),p,zm,'pchip');
    midz = sum(zm<=max(z_p(:,pidx)));
    for idz=midz+1:length(zm)
      p_zm(idz,pidx) = p_zm(idz-1,pidx)*(1-g*dz/2/RD/Tv_zm(idz-1,pidx))/(1+g*dz/2/RD/Tv_zm(idz-1,pidx));
    end
  end
  rho_env = p_z./Tv_env/RD;
  theta_e_env = eval_theta_e(T_env,q_env,0,0,p_z);
  rho_zm = p_zm./Tv_zm/RD;
  theta_e_zm = eval_theta_e(T_zm,q_zm,0,0,p_zm);
  theta_es_zm= eval_theta_es(T_zm,p_zm);
  theta_zm= eval_theta(T_zm,p_zm);
  rh_zm = q_zm./eval_qsat(T_zm,p_zm);
end

function [T_p, q_p, rho_p, z_p] = INTERP_P_GRID(plev, ta, hus)
global RD g p
  % Interpolate sounding values to 5-hPa p-grid (_p)
  T_p = nan([length(p),size(ta,2)]);
  q_p = nan([length(p),size(ta,2)]);
  Tv_p = nan([length(p),size(ta,2)]);
  rho_p = nan([length(p),size(ta,2)]);
  z_p = nan([length(p),size(ta,2)]);
  % pidx = find(prw==max(prw));
  for pidx=1:size(ta,2)
    T_p(:,pidx) = interp1(plev,ta(:,pidx),p,'pchip');
    q_p(:,pidx) = interp1(plev,hus(:,pidx),p,'pchip');
    % Fix negative q_p due to interp1
    q_p(q_p(:,pidx)<0,pidx) = 0;
    Tv_p(:,pidx) = eval_Tv(T_p(:,pidx),q_p(:,pidx));
    rho_p(:,pidx) = p./Tv_p(:,pidx)/RD;
    % Compute height assuming hydrostatic balance
    z_p(:,pidx) = (RD/g*(p(1)-p(2)))*cumtrapz(Tv_p(:,pidx)./p);
  end
end

function [plev, pr, ta, hus, date] = LOAD_SOUNDINGS(armbefn, arm)
  % Load sounding data
  if arm<=2 % Manus or Nauru
    plev = ncread(armbefn,'pressure'); % Pa
    pr = ncread(armbefn,'pr'); % mm/h
%     prw = ncread(armbefn,'prw'); % mm
    ta = ncread(armbefn,'ta'); % K
    hus = ncread(armbefn,'hus'); % kg/kg
%     rh = ncread(armbefn,'rh'); % in %; rh defined w.r.t. liquid (not ice)
    time = ncread(armbefn,'time_bnds'); % 'seconds since 1970-1-1 0:00:00 0:00'
    avail = sum(~isnan(ta),2)==37 & sum(~isnan(hus),2)==37 & ~isnan(pr);
    pr = pr(avail);
%     prw = prw(avail);
    ta = ta(avail,:)';
    hus = hus(avail,:)';
%     rh = rh(avail,:)';
    date = time(avail)/3600/24 + datenum(1970,1,1);
  elseif arm==3 % GOAmazon
    plev = ncread(armbefn,'pressure')*1e2; % hPa-->Pa
    pr = ncread(armbefn,'precip_rate_sfc'); % mm/h
    ta = ncread(armbefn,'temperature_p')'; % K
    rh = ncread(armbefn,'relative_humidity_p')'; % in %
    time = ncread(armbefn,'time'); % 'seconds since 2014-1-1 0:00:00 0:00'
    date = time/3600/24 + datenum(2014,1,1);
    hh = floor(mod(date,1)*24); 
    % GOAmazon hh==5 ---> Local Time 0030
    avail = sum(~isnan(ta),2)==37 & sum(~isnan(rh),2)==37 & ~isnan(pr) & hh==5;
    pr = pr(avail);
    ta = ta(avail,:)';
    rh = rh(avail,:)';
    hus = eval_q(eval_es(ta).*rh/100,repmat(plev,[1,size(ta,2)])); % Evaluated w.r.t. liquid water
    date = date(avail);
  end
end
