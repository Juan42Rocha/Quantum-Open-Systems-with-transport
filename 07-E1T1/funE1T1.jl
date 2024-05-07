include("fun.jl")


function SigZ(N)
  t = collect(N:-1:-N);
  return Diagonal(t)
end


function Jp(N,m)
    return sqrt(Complex(N*(N+1)-m*(m+1)))
end

function Jl(N,m)
    return sqrt(Complex(N*(N+1)-m*(m-1)))
end

function SigXY(N)
    m = collect(-N:1:N);
    t = Integer(2*N+1);
    im = collect(1:Integer(t));
    Jx = zeros(Complex,t,t);
    Jy = zeros(Complex,t,t); 
   
    for i in im
        for j in im
            jpij = 0;
            jlij = 0;
            if i==j+1
                jpij = Jp(N,im[i]);
            end
            if i==j-1
                jlij = Jl(N,im[i]);
            end
            Jx[i,j] = 0.5*(jpij+jlij); 
            Jy[i,j] = 0.5*(jpij-jlij);
        end
    end 
    return Jx,Jy
end
