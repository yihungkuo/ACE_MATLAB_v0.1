function wr = reconstruction_wr(k,v)
% wr defined by Eq. (2.58) in Shu (1998).
% Redundancy: length(v) = 2k-1.
  alphar = reconstruction_alphar(k,v);
  wr = alphar / sum(alphar);
end