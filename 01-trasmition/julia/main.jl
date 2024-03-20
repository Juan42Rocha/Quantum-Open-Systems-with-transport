include("fun.jl")

N     = 10

eps   = 0
gamma = 1
eta   = 1e-5
etas  = [1e-3,1]
Enum  = 500


function main(N,eps,gamma,etas,Enum)

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


for i in 1: etan;
  eta = etas[i];

  for Econ in 1:Enum;
    E = Eval[Econ];

    Sigs = eta*(E/2-1im*sqrt(1-(E/2)^2))*Nus;
    Sigd = eta*(E/2-1im*sqrt(1-(E/2)^2))*Nud;

    H = GeneraH(N,eps,gamma) + 1im*eta*(Nus+Nud) ; 
    H2 = -GeneraH(N,eps,gamma) ; 

    G = inv(E*Id-H); 
    G2 = GreenFun(E*Id,H2,Sigs+Sigd);#inv(E*Id - H2 - (Sigs+Sigd));

    A = -2*imag(G);
    
    A2 = -2*imag(G2);

    Rho[Econ,i] = tr(A2)/2/pi;


    TE[Econ,i]  =  real(4*tr((eta*Nus)*G*(eta*Nud)*G'));
    TE2[Econ,i] =  real(4*tr((imag(Sigs))*G2*(imag(Sigd))*G2'));

  end
  Rho[:,i] = Rho[:,i]./ findmax(Rho[:,i])[1];



end 

return TE2

end

