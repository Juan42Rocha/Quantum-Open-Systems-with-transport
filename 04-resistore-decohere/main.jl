##
#
# Esto es para introducir decoherencian el sistema
# la idea es que al con decoherencia la trasmisión del sistema
# se comporta como un sistema "real" y uno un conductor  perfecto
#
##

include("fun.jl")

function main()
Nr       = 500;  # relizaciones
Ndim     = 150;   # Max dimension
dcor     = 0.001;   # decore nu

step = 10


Nsize = collect(step:step:Ndim);

R   =zeros(Int(Ndim/step));

ind=1;
for i in Nsize;
    R[ind]= 1/Promedio_transmision_cerca_0(Nr,i,dcor);
    ind=ind+1;
    print(" ",ind*step)
end

io = open("nu"*string(dcor)*"-1.dat","w")
plot(Nsize,R)
#plot!(0:0.1:20, exp.(0.5collect(0:0.1:20)) )
#savedats(io,Nsize,R)


for i in 1:Int(Ndim/step)
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


