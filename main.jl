using LinearAlgebra
using DelimitedFiles



function main();

  g = 1 ;
  eA = zeros(Complex,4,4);
  l1 = -2*g; 
  l2 = l1 + 2*im;
  l3 = l1 - 2*im;
 
  st0 =  [0 ; 0 ; 1];
  st = zeros(Complex,3,101);
 
  
  for it in collect(1:101);
    T  = collect(0:10/100:10); 
    t = T[it];

    t1 = -exp(l1*t) + exp(l2*t) -2*exp(l2*t);
    t2 = 2*im*(exp(l2*t) - exp(l3*t));
    t3 = -exp(l1*t) - exp(l2*t) +2*exp(l2*t);   
    t4 = -4*exp(l2*t) + 2*exp(l3*t);
 
    eA = (-1/2).*[t1  t2  t3 ;
                 t2  t4  -t2;
                 t3 -t2  t1];
    @show t
    IeAb = -(g*l1).*[1 ; 0 ; 1] .* (1-exp(l1*t)); 
    st[:,it] = complex(eA)*complex(st0) + complex.(IeAb);     

  end
  
  open("testRe.dat","w") do io
    writedlm(io,real.(st)')
  end

  open("testIm.dat","w") do io
    writedlm(io,imag.(st)')
  end

end
