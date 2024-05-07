##
#
# Esto es para introducir decoherencian el sistema
# la idea es que al con decoherencia la trasmisión del sistema
# se comporta como un sistema "real" y uno un conductor  perfecto
#
##

include("fun.jl")

function main()
Nr       = 5000;  # relizaciones
Ninit    = 2 ;
Ndim     = 30;   # Max dimension
dcor     = 0.5;   # decore nu

step = 2


Nsize = collect(Ninit:step:Ndim);

R   =zeros(Int((Ndim-Ninit)/step)+1);

ind=1;
for i in Nsize;
    R[ind]= 1/Promedio_transmision_cerca_0(Nr,i,dcor);
    ind=ind+1;
    print(" ",(ind-1)*step+Ninit)
end

io = open("Cnu"*string(dcor)*".dat","w")
#plot(Nsize,R)
#plot!(0:0.1:20, exp.(0.5collect(0:0.1:20)) )
#savedats(io,Nsize,R)


for i in 1:Int((Ndim-Ninit)/step + 1)
  @printf(io,"%f   %f \n",Nsize[i],R[i])
end

close(io)


end











#
# io = open("filename.dat","w")
#  savedats
#
#  close(io)

function savedats(io,x,y)

  n = size(x)
  for i=1:n
     @printf(io,"%f   %f ",x[i],y[i])
  end
  @printf(io,"\n")
  

end


