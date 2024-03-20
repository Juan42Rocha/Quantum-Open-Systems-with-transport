using LinearAlgebra
using Random


function buildH(Ndim)
    Hmat = zeros(Float64, (Ndim, Ndim))
    for i in 1:1:Ndim-1
        Hmat[i,i+1] = 1.0
        Hmat[i+1,i] = 1.0
        Hmat[i,i] = randn()
    end 
    Hmat[Ndim, Ndim] = randn()
    return Hmat
end

function buildSigma(Ndim, npos, eta)
    Sig = zeros(Complex{Float64}, (Ndim,Ndim))
    Sig[npos,npos] = 1im*eta
    return Sig
end 

function Tvec1(Ndim, SigR, Sigj, contacts, G)
    Tarr = zeros(Float64, size(contacts)[1])
    tt = 1
    for ii in contacts
        Tarr[tt] = 4.0*real(tr(imag(SigR)*G*imag(Sigj[tt])*G'))
        tt+=1
    end 
    return Tarr
end

function Tvec2(Ndim, SigR, Sigj, contacts, G)
    Tarr = zeros(Float64, size(contacts)[1])
    tt = 1
    for ii in contacts
        Tarr[tt] = 4.0*real(tr(imag(Sigj[tt])*G*imag(SigR)*G'))
        tt+=1
    end	
    return Tarr
end

function buildtau(Ndim, Sigj, SigS, SigD, contacts, G)
    tau = zeros(Float64, (size(contacts)[1], size(contacts)[1]))
    ss = 1
    for ii in contacts
        for jj in contacts
            tt = 1
            if ss!=tt
                tau[ss,tt] = -4.0*real(tr(imag(Sigj[ss])*G*imag(Sigj[tt])*G'))
            end
            tt+=1
        end
        Tsum = 0.0
        Nsize = size(contacts)[1]
        for kk in 1:1:(Nsize+2)
            if kk==(Nsize+1)
                Tsum += 4.0*real(tr(imag(Sigj[ss])*G*imag(SigS)*G'))
            elseif kk==(Nsize+2)
                Tsum += 4.0*real(tr(imag(Sigj[ss])*G*imag(SigD)*G'))
            elseif kk==ss
                continue
            else
                Tsum += 4.0*real(tr(imag(Sigj[ss])*G*imag(Sigj[kk])*G'))
            end
        end
        tau[ss,ss] = Tsum
        ss+=1
    end
    return tau 
end

function main()
    eta = 1.0; nu = 0.5; pnu = 0.7
    E = 0.00001
    numrealiz = 200
    Nmin = 10; Nmax = 30
    ss = 1
    printf = open("rescenband_Nmin$(Nmin)Nmax$(Nmax)numrealiz$(numrealiz)eta$(round(eta, digits=2))nu$(round(nu,digits=2))pnu$(round(pnu,digits=2)).dat", "w")
    for Ndim in Nmin:2:Nmax##[10,20,30,40,50,60,70] ##Nmin:5:Nmax
        Transm = 0.0
        SigS = buildSigma(Ndim, 1, eta)
        SigD = buildSigma(Ndim, Ndim, eta)
        # Build all the matrices for the ficticious sites
        Sigj = Matrix{Complex{Float64}}[]
        contacts = Array{Int64}[]
        for site in 1:1:Ndim
           if rand() <= pnu
               push!(Sigj, buildSigma(Ndim, site, nu))
               push!(contacts, [site])
           end 
        end 
        @show Ndim, contacts
        # Run over all the realizations
        for irealiz in 1:1:numrealiz   
            Hmat = buildH(Ndim)
            Hmat = Hmat + SigS + SigD
            for site in contacts
                sitei = site[1]
                Hmat[sitei,sitei] += 1im*nu
            end
            G = inv(E*I - Hmat)
            # Compute transmission vectors
            TvecD = Tvec1(Ndim, SigD, Sigj, contacts, G)
            TvecS = Tvec2(Ndim, SigS, Sigj, contacts, G)
            # Compute matrix tau
            tau = buildtau(Ndim, Sigj, SigS, SigD, contacts, G)
            Transm += 4.0*real(tr(imag(SigD)*G*imag(SigS)*G')) + TvecD'*inv(tau)*TvecS
        end
        Transm = Transm/numrealiz
        println(printf, Ndim, " ", 1.0/Transm)
        ss+=1 
    end
    close(printf)
end # main

main()
