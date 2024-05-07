include("funE1T1.jl")

 # Se busca calcula la tranmision entre con el 
 # hamiltonia no con la forma 
 # 
 #  H =(eps_0 - i gamma/2)I +2k Sx - i gamma Sz
 #
 # donde I es la matriz identidad 
 # Sx es obien  el operador de mometno angular en x o un escalar por la matiz de pauli 
 # Sz lo mismo pero en z 
 # aqui se toman como 
 #   Sx(i,j) =  delta(i-1,j)+delta(i+1,j)
 #   sz(i,j) =  delta(i,j)m_i



function main()
    ep = 0;  gamma = 0.1; k = 0.5;  N = 3;

    Hr = GeneraH(N,ep,2*k/sqrt(2));
    Hi = GeneraH(N,-im*0.5*gamma,0);
    H1 = -im*gamma*SigZ((N-1)/2) ;
    Heff = Hr+Hi+H1;
    
    Eival, Eivec = eigen(Heff);    
 
    eio  = open("eigN"*string(N)*"g"*string(gamma)*"k"*string(k)*".dat","w")
    for k in 1:N
        @printf(eio,"%f  %f  \n",real(Eival[k]),imag(Eival[k]) )
    end 
    close(eio)

    nE = 1000
   

    Evals = collect(-2.1:4.1/nE:2);
  #  T = zeros(nE+1);  
    io = open("2T_N"*string(N)*"g"*string(gamma)*"k"*string(k)*".dat","w")
    for Ei  in 1:nE+1
        E = Evals[Ei];
        G = GreenFun(E*I,Hr,H1+Hi);

        TbS = Mat_TransmisionS(Heff[1,1],Heff[2,2],G,1,N);
        TDa = Mat_TransmisionS(Heff[N,N],Heff[2,2],G,N,N);
        tau = Mat_Tau(Heff[1,1],Heff[N,N],Heff[2,2],G,N);
        println(det(tau)," ",E," ",Ei)
        T = Calcula_transmision(Heff[1,1],Heff[N,N],G,1,N)
        
invtau = inv(tau);
        temp =  TDa'*invtau*TbS;
        
        T2 = T + temp[1];
        @printf(io,"%f  %f  %f \n",E,T,T2)

       
    end 
    println()
    close(io)   
return (Heff)

end
a = main()
