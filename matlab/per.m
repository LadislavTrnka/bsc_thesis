% Permutation matrix
% To generate a permutation matrix, which changes rows and columns of the Chebyshev grid in a vector of the considered discrete function

function P = per(NXN,NYN)
    P = zeros(NXN*NYN);
    for m = 0:NXN-1
      for k = 0:NYN-1
           P(k+1+NYN*m,NXN*k+1+m)=1;
      end  
    end
end

