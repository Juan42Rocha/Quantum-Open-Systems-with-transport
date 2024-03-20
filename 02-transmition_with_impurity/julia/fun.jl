using LinearAlgebra
using Plots

function  GeneraH(N,eps,gamma)

  H = zeros(N,N)


  for i in 2:N;
    H[i,i-1] = gamma;
    H[i,i]   = eps/2;
  end
  H = H + H' 
  H[1,1] = eps
  return H
end


function AddImpurtiy(H,n,t)
  N = size(H,2);
  newH = zeros(N+1,N+1);
  newH[1:end-1,1:end-1] = H;
  newH[end,n] = t;
  newH[n,end] = t;


return newH
end

function expan1(M)
  N = size(M,2);
  nM = zeros(Complex,N+1,N+1);
  nM[1:end-1,1:end-1] = M;


return nM
end


function GreenFun(E,H,SIGMA)

  G = inv(E-H-SIGMA);

  return G

end
