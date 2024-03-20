using LinearAlgebra
using Random
using Distributions

include("../../fun.jl")

function buildH(Ndim)
    Hmat = zeros(Float64, (Ndim, Ndim))
    dummy = Normal(0.0,2.0)
    for i in 1:1:Ndim-1
        Hmat[i,i+1] = 1.0
        Hmat[i+1,i] = 1.0
        Hmat[i,i] = rand(dummy)
    end 
    Hmat[Ndim, Ndim] = rand(dummy)
    return Hmat
end

function buildSigma(Ndim, npos, eta)
    Sig = zeros(Complex{Float64}, (Ndim,Ndim))
    Sig[npos,npos] = 1im*eta
    return Sig
end 

function Tvec1(Ndim, SigR, Sigj, G)
    Tarr = zeros(Float64, Ndim)
    for ii in 1:1:Ndim
        #Tarr[ii] = 4.0*real(tr(imag(SigR)*G*imag(Sigj[ii])*G'))
        Tarr[ii] = 4.0*real(tr(imag(SigR)*G*imag(Sigj[ii])*G'))
    end 
    return Tarr
end

function Tvec2(Ndim, SigR, Sigj, G)
    Tarr = zeros(Float64, Ndim)
    for ii in 1:1:Ndim
        Tarr[ii] = 4.0*real(tr(imag(Sigj[ii])*G*imag(SigR)*G'))
    end	
    return Tarr
end

function buildtau(Ndim, Sigj, SigS, SigD, G)
    tau = zeros(Float64, (Ndim, Ndim))
    for ii in 1:1:Ndim
        for jj in 1:1:Ndim
            if ii!=jj
                tau[ii,jj] = -4.0*real(tr(imag(Sigj[ii])*G*imag(Sigj[jj])*G'))
            end
        end
        Tsum = 0.0
        for kk in 1:1:(Ndim+2)
            if kk==(Ndim+1)
                Tsum += 4.0*real(tr(imag(Sigj[ii])*G*imag(SigS)*G'))
            elseif kk==(Ndim+2)
                Tsum += 4.0*real(tr(imag(Sigj[ii])*G*imag(SigD)*G'))
            elseif kk==ii
                continue
            else
                Tsum += 4.0*real(tr(imag(Sigj[ii])*G*imag(Sigj[kk])*G'))
            end
         end
         tau[ii,ii] = Tsum
    end
    return tau 
end

function main()
    eta = 1.0; nu = 0.2
    E = range(-2.0,2.0,200)
    numrealiz = 50
    Nmin = 10; Nmax = 30
    printf = open("rescenband_Nmin$(Nmin)Nmax$(Nmax)numrealiz$(numrealiz)eta$(round(eta, digits=2))nu$(round(nu,digits=2))sigma2.0.dat", "w")
    for Ndim in [15]
        Transm = zeros(Float64, 200)
        SigS = buildSigma(Ndim, 1, eta)
        SigD = buildSigma(Ndim, Ndim, eta)
        # Build all the matrices for the ficticious sites
        Sigj = Matrix{Complex{Float64}}[]
        for site in 1:1:Ndim
           push!(Sigj, buildSigma(Ndim, site, nu))
        end 
        # Run over all the realizations
        for irealiz in 1:1:numrealiz   
            Hmat = buildH(Ndim)
            Hmat = Hmat + SigS + SigD
            for site in 1:1:Ndim
                Hmat[site,site] += 1im*nu
            end
            tt=1
            for Ei in E
                @show Ndim, irealiz, tt
                G = inv(Ei*I - Hmat)
                # Compute transmission vectors
                TvecD =  Mat_TransmisionD(im*eta,im*nu,G,1,Ndim);
                #Tvec1(Ndim, SigD, Sigj, G)
                TvecS = Mat_TransmisionS(im*eta,im*nu,G,Ndim,Ndim) # Tvec2(Ndim, SigS, Sigj, G)
                # Compute matrix tau
                tau = Mat_Tau(im*eta,im*nu,G,Ndim); #buildtau(Ndim, Sigj, SigS, SigD, G)
                invtau = inv(tau);
                temp = TvecD'*invtau*TvecS
                
                Transm[tt] += Calcula_transmision(im*eta,im*eta,G,1,Ndim)+temp[1,1];#4.0*real(tr(imag(SigD)*G*imag(SigS)*G')) + TvecD'*inv(tau)*TvecS
                tt+=1
            end
        end
        Transm = Transm/numrealiz
        ss=1
        for Ei in E
            println(printf, Ei, " ", Transm[ss])
        ss+=1
        end
    end
    close(printf)
end # main

main()
