function crj = reconstruction_cij(k,r)
% crj defined by Eq. (2.21) in Shu (1998):
%   V(i+1/2) ≈ ∑_{j=0}^{k-1} crj V(i-r+j),
% where k>0 is oder of accuracy, r non-negative integer left shift.
  crj = zeros(k,1);
  for j=0:k-1
    for m=j+1:k
    	deno = 1.0;
      for l=0:k
        if l~=m
          deno = deno * (m-l);
        end
      end
      nume = 0.0;
      for l=0:k
        if l~=m
          temp = 1.0;
          for q=0:k
            if q~=m && q~=l
              temp = temp * (r-q+1);
            end
          end
          nume = nume + temp;
        end
      end
      crj(j+1) = crj(j+1) + nume/deno;
    end
  end
end