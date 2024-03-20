##
#
# Esto es para la locación de anderso
#
##

include("fun.jl")

function main()
Nr=500; #relizaciones
Ndim=100; #dimension

R=zeros(20);
ind=1;
for i in 5:5:Ndim 
    R[ind]= 1/Promedio_transmision_cerca_0(Nr,i);
    ind=ind+1;
    print(ind)
end



plot(1:1:20,R,yaxis=:log)
plot!(0:0.1:20, exp.(0.5collect(0:0.1:20)) )



end
