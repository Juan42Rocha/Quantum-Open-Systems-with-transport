include("fun.jl")
include("functions.jl")


function main(gamma,Nstep,Rel,np)


st0     = [0;1];
#gamma   = 1;
#Nstep = Int(1e3);
#Rel     =  100;
dt      = np*pi/(Nstep+1.0) #1e-5*2*pi;

sx = [0 1 ; 1 0];
c = sqrt(gamma)*[0 1 ; 0 0];


H = sx - 1im/2 *c'*c;
println(H);

#st = OneRel(st0,gamma,dtH,c,Nstep);
rho = MultiRel(st0,gamma,dt,H,c,Nstep,Rel);


A = [ 0   im  -im  gamma ;
      im  -0.5*gamma  0  -im;
      -im  0  -0.5*gamma  im;
      0   -im  im    -gamma];


rho2 = rk4(A,reshape(st0*st0',4,1),Nstep,dt);


#plot(1:Nstep,abs.(st[1,:]).^2)
t = collect(0:Nstep)*np/(Nstep+1)
d = sqrt(64.0-gamma^2); 
#rho3 = 0.5 .+0.5*(abs(st0[1,1])^2 -abs(st0[2,1])^2).*cos.(2*t*pi) #abs(st0[2,1])^2 .*exp.(-0.75*gamma.*t).*cos.(d.*t./4);
rho3 = (abs(st0[2,1])^2 .*exp.(-0.75*gamma.*t*pi).*cos.(d.*t*pi./4))*0.5 .+0.5;
return t,rho,rho2,rho3;

end

