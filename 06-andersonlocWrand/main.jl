include("fun.jl")

N     = 10

eps   = 0
gamma = 1
eta   = 1e-5
etas  = [1]
Enum  = 500


function main(N,Enum,dcor)
  etas = [0.999999];
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



  for Econ in 1:Enum;
    E = Eval[Econ];
    #Si quieres distribucion normal
    d=Normal(0,);
    eta=0.9999*((E/2) - im*sqrt(1 - ((E/2)^2)));

    eps=Diagonal(rand(d,N));
    #eps=Diagonal(rand(N));

    Id=eye(N);
    Sigs = zeros(N,N);
    Sigd = zeros(N,N);
    Sigs[1,1] = 1;
    Sigd[N,N] = 1;
    Mtau = zeros(N,N);
    

    #G = GreenFun(E*Id,H2,eps+Sigs+Sigd);#inv(E*Id - H2 - (Sigs+Sigd));
    G = GreenFun(E*Id, GeneraH(N,0,1)+Sigmadecore(N,dcor),eps+eta*(Sigs+Sigd));
  
    if dcor != 0 ;
      TbS =  Mat_TransmisionS(eta,1im*dcor,G,1,N);
      TDa =  Mat_TransmisionS(eta,1im*dcor,G,N,N);
      Mtau = Mat_Tau(eta,1im*dcor,G,N);
      invtau = inv(Mtau);
      temp = TDa'*invtau*TbS
      TE[Econ,i] = Calcula_transmision(eta,eta,G,1,N) + temp[1,1]; 
    else
      TE[Econ,i] Calcula_transmision(eta,eta,G,1,N) # + temp[1,1]; 
    end

end 

return TE2

end

