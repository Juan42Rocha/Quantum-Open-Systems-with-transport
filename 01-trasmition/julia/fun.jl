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


function GreenFun(E,H,SIGMA)

  G = inv(E-H2-SIGMA);

  return G

end
