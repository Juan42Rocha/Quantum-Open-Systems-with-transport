exec("test.sci",-1)

N     = 10

eps   = 0
gamma = 1
eta   = 1e-5
etas  = [1e-3]
etan  = length(etas);
Enum  = 500;
Eval  = [-2:4/Enum:2]


Nus = zeros(N,N);
Nud = zeros(N,N);
Sigs = zeros(N,N);
Sigd = zeros(N,N);
Nus(1,1) = 1;
Nud(N,N) = 1;

Id = eye(N,N);
Rho = zeros(Enum,etan);
TE = zeros(Enum,etan);
TE2 = zeros(Enum,etan);


for i = 1: etan
  eta = etas(i)

  for Econ = 1:Enum
    E = Eval(Econ);

    Sigs = (E/2-%i*sqrt(1-(E/2)^2))*Nus;
    Sigd = (E/2-%i*sqrt(1-(E/2)^2))*Nud;

    H = GeneraH(N,eps,gamma) - %i*eta*(Nus+Nud) ; 
    H2 = GeneraH(N,eps,gamma) ; 

    G = inv(E*Id-H); 
    Gp = (E*Id -H2 -%i*eta*(Sigs+Sigd));
    G2 = double( invr(E*Id-H2-%i*eta*(Sigs+Sigd)));


    A = -2*imag(G);
    
    A2 = -2*imag(G2);

    Rho(Econ,i) = trace(A2)/2/%pi;
   // if (Rho(Econ,i)<0)
   // pause
   // end


    TE(Econ,i)  =  real(4*trace((eta*Nus)*G*(eta*Nud)*G'));
    TE2(Econ,i) =  real(4*trace((imag(eta*Sigs))*G2*(imag(eta*Sigd))*G2'));

  end
  Rho(:,i) = Rho(:,i)./ max(Rho(:,i));

 // [eigVec,t1]=spec(H);
 // eigVal = diag(t1);q

end 
disp(TE2)


//exitq
