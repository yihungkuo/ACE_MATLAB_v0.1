function v_weno = reconstruction_weno(k,v)
% Reconstruction using Eq. (2.52) in Shu (1998):
% 
% 	V(i+1/2) ≈ ∑_{r=0}^{k-1} wr Vr(i+1/2),
% 
% Redundancy: length(v) = 2k-1 is the order of accuracy.
% Note that v = V(i-k+1:i+k-1).
  wr = reconstruction_wr(k,v);
  v_weno = 0;
  for r=0:k-1
    v_weno = v_weno + wr(r+1) * dot(reconstruction_cij(k,r),v(k-r:end-r));
  end
end