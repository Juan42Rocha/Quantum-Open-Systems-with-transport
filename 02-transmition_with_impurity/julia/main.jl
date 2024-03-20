include("fun.jl")

N     = 10

eps   = 0
gamma = 1
eta   = 1e-5
etas  = [1e-3,1]
Enum  = 500


function main(N,eps,gamma,etas,Enum, nimpu, nt)

etan  = length(etas);
Eval  = LinRange(-2,2,Enum);


Nus = zeros(N,N);
Nud = zeros(N,N);
Sigs = zeros(N,N);
Sigd = zeros(N,N);
Nus[1,1] = 1;
Nud[N,N] = 1;

Id = diagm(ones(N));
Rho = zeros(Enum,etan);
TE = zeros(Enum,etan);
TE2 = zeros(Enum,etan);
TE3 = zeros(Enum,etan);
G3 =  zeros(N,N);

for i in 1: etan;
  eta = etas[i];

  for Econ in 1:Enum;
    E = Eval[Econ];

    Sigs = eta*(E/2-1im*sqrt(1-(E/2)^2))*Nus;
    Sigd = eta*(E/2-1im*sqrt(1-(E/2)^2))*Nud;

    nSigs = expan1(real(Sigs));
    nSigd = expan1(real(Sigd));
    nId = expan1(Id);

    H = GeneraH(N,eps,gamma) + 1im*eta*(Nus+Nud) ; 
    H2 = -GeneraH(N,eps,gamma) ; 
    H3 = AddImpurtiy(H2,nimpu,nt);

 #   G = inv(E*Id-H); 
 #   G2 = GreenFun(E*Id,H2,Sigs+Sigd);#inv(E*Id - H2 - (Sigs+Sigd));
    G3 = GreenFun(E*nId,H3,nSigs+nSigd);#inv(E*Id - H2 - (Sigs+Sigd));

  #  A2 = -2*imag(G2);
 #   Rho[Econ,i] = tr(A2)/2/pi;

 #   TE[Econ,i]  =  real(4*tr((eta*Nus)*G*(eta*Nud)*G'));
 #   TE2[Econ,i] =  real(4*tr((imag(Sigs))*G2*(imag(Sigd))*G2')); 
# println(imag(nSigs)*G3*imag(nSigd)*G3' );
# println(real( 4*tr(imag(nSigs)*G3*imag(nSigd)*G3' )) );
 # TE3[Econ,i] = real( 4*tr(imag(nSigs)*G3*imag(nSigd)*G3' ));
 TE3[Econ,i]  =  4*imag(nSigs+nsigd)*G3[1,end-1];
  end
  #Rho[:,i] = Rho[:,i]./ findmax(Rho[:,i])[1];



end 

return G3

end

