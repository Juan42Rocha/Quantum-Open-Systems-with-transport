# ========================================
# funciones para intentar resolver e. maestras
# (d/dt) vec[rho] =  A vec[rho]
# (d/dt) vec[rho] =  f( vec[rho] )
# ========================================


function  rk4step(A,rho0,dt)

  k1 = A*rho0;

  k2 = A*(rho0 + 0.5*dt*k1);
  
  k3 = A*(rho0 + 0.5*dt*k2);

  k4 = A*(rho0 + dt*k3);

  rho1 = rho0 + dt*(k1 + 2*k2 + 2*k3 + k4)/6;                              

  return rho1;
end

 
function rk4(A,rho0,N,dt)

  rho = zeros(Complex,size(rho0,1),N+1);
  rho[:,1] = rho0  ;
 
  for k =1:N
    rho[:,k+1] = rk4step(A,rho[:,k],dt);
  end

  return rho;
end




