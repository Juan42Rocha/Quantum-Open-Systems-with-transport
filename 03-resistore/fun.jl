using LinearAlgebra
using Random
using Plots
using Distributions

#
# Crea matriz tridiagonal 
#
#
function  GeneraH(N,eps,gamma)

  H = zeros(N,N)

  for i in 2:N;
    H[i,i-1] = gamma;
    H[i,i]   = eps/2;
  end
  H = H + H' ;
  H[1,1] = eps;
  return H
end

#funcion dos para decoerencias
function  Sigmadecore(N,nu)

  H = zeros(Complex,N,N);

  for i in 1:N;
    H[i,i]   =  im *nu;
  end
  return H
end

function GreenFun(E,H,SIGMA)
  G = inv(E-H-SIGMA);
  return G
end

function eye(n)
  return Matrix(I,n,n)
end

function transmision_etavar_rand(E,N)

    #Si quieres distribucion normal
    #d=Normal();
    

    #eps=Diagonal(rand(d,N));
    eps=Diagonal(rand(N));

    eta=0.9999*((E/2) - im*sqrt(1 - ((E/2)^2)));

    Id=eye(N);
    Sigs = zeros(N,N);
    Sigd = zeros(N,N);
    Sigs[1,1] = 1;
    Sigd[N,N] = 1;
    
    G = GreenFun(E*Id, GeneraH(N,0,1)+Sigmadecore(N,0.5),eps+eta*(Sigs+Sigd));
    
    return real(4*tr(imag.(eta*Sigs)*G*imag.(eta*Sigd)* G' ) )
end


function Promedio_transmision_cerca_0(Nr,N)
    T=0.0;
    
    for i in 1:1:Nr
        T=T+ transmision_etavar_rand(0.0,N)
    end
    
   return T/Nr
end
