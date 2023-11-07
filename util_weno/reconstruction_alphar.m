function alphar = reconstruction_alphar(k,v,epsilon)
% alphar defined by Eq. (2.59) in Shu (1998).
% Redundancy: length(v) = 2k-1.
  if nargin==2
    epsilon = 1e-10;
  end
  alphar = reconstruction_dr(k)./ (epsilon + reconstruction_betar(k,v)).^2;
end