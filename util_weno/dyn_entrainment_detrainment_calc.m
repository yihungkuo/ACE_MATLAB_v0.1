function dMFdz = dyn_entrainment_detrainment_calc(MF)
global k dz 
  [Nzm, Nt] = size(MF);
  dMFdz = nan(Nzm,Nt);

  MFext = [-MF(k:-1:2,:); MF];
  
  Nzext = Nzm+k-1; 
  dMFextdz = zeros(Nzext,1);
  LMFh = zeros(Nzm,1); % MFh from more left stencil points
  RMFh = zeros(Nzm,1); % MFh from more right stencil points
  for idt=1:Nt
    % Compute MFh(zm+dz/2) from MF(zm) as the mean of up- and down-wind approximations 
    for idz=k+1:Nzext-k+1
      LMFh(idz-k) = reconstruction_weno(k, MFext([idz-k:idz+k-2],idt) );
      RMFh(idz-k) = reconstruction_weno(k, MFext(flip([idz-k+1:idz+k-1]),idt) );
    end
    MFh = (LMFh + RMFh)/2;
    for idz=k+1:Nzext-k
      % d(MF=rho_0*w)/dz
      dMFextdz(idz)...
          = ( reconstruction_weno(k, MFext( ([idz-k+1:idz+k-1])*(MFh(idz-k+1)>=0)+flip([idz-k+2:idz+k])*(MFh(idz-k+1)<0) ,idt) ) ...
             - reconstruction_weno(k, MFext( ([idz-k:idz+k-2])*(MFh(idz-k)>=0)+flip([idz-k+1:idz+k-1])*(MFh(idz-k)<0) ,idt) ) )/dz;
    end%idz
    % First zero at zm=0 (dXdz~=0 but dXdt=0 at the end of computation)
    % zeros(k,1) for domain top<---to be improved!!!!!
    dMFdz(:,idt) = [0; dMFextdz(k+1:Nzext-k); zeros(k,1)]; % Positive/negative=entrainment/detrainment
  end%idt
end