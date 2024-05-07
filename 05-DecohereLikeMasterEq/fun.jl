using Random
using LinearAlgebra
#using Plots
using Printf

function StepJump(H,st,gamma,dt,c);
  dp = gamma*abs(st[2,1])^2*dt;
  nst = zeros(Complex,2,1);
  if  rand() < dp  ;
   # println("jump")
    nst = c*st;
    nst = nst./norm(nst);
  else
    nst = st;
    nst =(I(2)-1im*dt*H)*st;
  end

  return nst;
end


function OneRel(st,gamma,dt,H,c,N);
 # (nst,temp) = size(st);
  stT = zeros(Complex,2,N+1);
  stT[:,1] = st;
#  println(st,stT[:,1]);
  for k in 1:N
   # stT[:,k+1] =(I(2)-1im*dt*H)*stT[:,k];
    stT[:,k+1] = StepJump(H,stT[:,k],gamma,dt,c);
    stT[:,k+1] =  stT[:,k+1]./norm(stT[:,k+1]);
  end

  return stT;
end


function MultiRel(st,gamma,dt,H,c,N,rel);

  rho = zeros(Complex,4,N+1);
  stR = zeros(Complex,2,N+1);
  for i in 1:rel
    stR =  OneRel(st,gamma,dt,H,c,N);
    rho = rho + makeRho(stR,N);

  end 

  return rho./rel;


end


function makeRho(st,N);
  rho =  zeros(Complex,4,N+1);
  for k in 1:N+1;
    rho[:,k] = reshape(st[:,k]*st[:,k]',4,1);
  end
  return rho;

end



function CM(A,B)
  return A*B-B*A
end
1
function aCM(A,B)
  return A*B+B*A
end

