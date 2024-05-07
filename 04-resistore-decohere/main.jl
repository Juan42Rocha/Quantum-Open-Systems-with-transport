##
#
# Esto es para introducir decoherencian el sistema
# la idea es que al con decoherencia la trasmisión del sistema
# se comporta como un sistema "real" y uno un conductor  perfecto
#
##

include("fun.jl")

function main()
<<<<<<< HEAD
Nr       = 5000;  # relizaciones
Ninit    = 2 ;
Ndim     = 30;   # Max dimension
dcor     = 0.5;   # decore nu

step = 2
=======
Nr       = 200;  # relizaciones
Ndim     = 30;   # Max dimension
Ninit    = 10;
dcor     = 0.5;   # decore nu

step = 1;
>>>>>>> e573dd039a6494128a48da8691763c4474b5b35b


Nsize = collect(Ninit:step:Ndim);

<<<<<<< HEAD
R   =zeros(Int((Ndim-Ninit)/step)+1);
=======
R   =zeros(size(Nsize,1));
>>>>>>> e573dd039a6494128a48da8691763c4474b5b35b

ind=1;
for i in Nsize;
    R[ind]= 1/Promedio_transmision_cerca_0(Nr,i,dcor);
    ind=ind+1;
<<<<<<< HEAD
    print(" ",(ind-1)*step+Ninit)
end

io = open("Cnu"*string(dcor)*".dat","w")
#plot(Nsize,R)
=======
    print(" ",ind*step+8)
end

io = open("Repli-nu"*string(dcor)*"-1.dat","w")
plot(Nsize,R)
>>>>>>> e573dd039a6494128a48da8691763c4474b5b35b
#plot!(0:0.1:20, exp.(0.5collect(0:0.1:20)) )
#savedats(io,Nsize,R)


<<<<<<< HEAD
for i in 1:Int((Ndim-Ninit)/step + 1)
=======
for i in 1:Int((Ndim-Ninit)/step)
>>>>>>> e573dd039a6494128a48da8691763c4474b5b35b
  @printf(io,"%f   %f \n",Nsize[i],R[i])
end

close(io)


end











#
# io = open("filename.dat","w")
#  savedats
#
#  close(io)

#function savedats(io,x,y)
#
#  n = size(x)
#  for i=1:n
#     @printf(io,"%f   %f ",x[i],y[i])
#  end
#  @printf(io,"\n")
#  
#
#end


