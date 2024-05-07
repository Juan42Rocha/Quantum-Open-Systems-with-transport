include("fun.jl")
include("functions.jl")


function main(st0,om,gamma,Nstep,Rel,np)


#st0     = [0;1];
#gamma   = 1;
#Nstep = Int(1e3);
#Rel     =  100;
dt      = np*pi/(Nstep+1.0); #1e-5*2*pi;

sx = [0 1 ; 1 0];
c = sqrt(gamma)*[0 1 ; 0 0];


H = sx.*om - 1im/2 *c'*c;
println(H);

#st = OneRel(st0,gamma,dtH,c,Nstep);
rho = MultiRel(st0,gamma,dt,H,c,Nstep,Rel);


A = [ 0   im*om  -im*om  gamma ;
      im*om  -0.5*gamma  0  -im*om;
      -im*om  0  -0.5*gamma  im*om;
      0   -im*om  im*om    -gamma];


rho2 = rk4(A,reshape(st0*st0',4,1),Nstep,dt);



#plot(1:Nstep,abs.(st[1,:]).^2)
t = collect(0:dt:pi*np);
 

rho3 = Analitica(t,st0*st0',om,gamma);  #(abs(st0[2,1])^2 .*exp.(-0.75*gamma.*t*pi).*cos.(d.*t*pi./4))*0.5 .+0.5;

file = open("test.dat","w")
for k in 1:Nstep+1
   @printf(file,"%f    %f    %f     %f\n",t[k],real(rho[4,k]),real(rho2[4,k]),real(rho3[k][1]))

end 

close(file)
return t,rho,rho2,rho3;

end


function Analitica(T,rho0,om,gamma)

    p0 = [rho0[2,2] ; imag(rho0[2,1])];
    
    lp = -3*gamma/4 + sqrt(complex(gamma^2-64*om^2))/4;
    lm = -3*gamma/4 - sqrt(complex(gamma^2-64*om^2))/4;
    
    gp = (gamma + sqrt(complex(gamma^2 -64*om^2)))/om;
    gm = (gamma - sqrt(complex(gamma^2 -64*om^2)))/om;
    
    c = om/(sqrt(complex(gamma^2-64*om^2))*2);
    

    p = [];
    for t in T
        M = [ -(gm*exp(lp*t)-gp*exp(lm*t))   8*(exp(lp*t)-exp(lm*t));
             -16*(exp(lp*t)-exp(lm*t))       gp*exp(lp*t)-gm*exp(lm*t)  ];
       # @show M,c, gp-gm, 1/c (gp-gm)*c
        B1 = [8 ; gm].*om*(1-exp(lm*t))./lm;
        B2 = [8 ; gp].*om*(1-exp(lp*t))./lp;
        push!(p,transpose( c.*(M*p0+B1-B2 ) ));
    end

    return p;          
end

st0   = [0;1];
om    = 5;
gamma = 0.4;
np    = 7;
Nstep = 2000;
Rel   = 2000;

t,r,r2,r3 = main(st0,om,gamma,Nstep,Rel,np);
@show st0
