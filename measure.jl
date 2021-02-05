function measure(sector, L,d,J,B,Tmin,Tmax,Tstep,runs)

    # average magnetization, average energy and extracting Cv and susceptibility as a function of T
    # Cv=(mean(E^2)-(mean(E))^2)/T^2
    # Xi=L*L*(mean(M.^2)-(mean(abs.(M)))^2)/T 
    # actual susceptibility requires to multiply with number of spins!
    # notice abs(M) - essential for correct behaviour of Xi for scaling, due to SSM happening for real
    # in small lattices
    # sector: indicates which Wilson loop sector is addresssed (+-1 across d dimensions)

    Cutoff=10000


    Temperature = Float32[]
    SpecificHeat=Float32[]
    MagnetizationSquared = Float32[]
    EnergyMean = Float32[]
    Susceptibility = Float32[]
    Wilson_loops_mean=Array{Int32,2}(undef,0,d) # average of array of wilson loops in directions 1,2,...d at each temperature step
    wilson_mean=ones(Float32,d)
    
    no_steps = Integer(1+round((Tmax-Tmin+Tstep)/Tstep))        
    p=Progress(no_steps,1.0)
    
    sigma=setlattice(L,d,sector);
    

    for T=Tmin:Tstep:Tmax 
        energy, magnetization,sigma, wilson_loops=IsingGaugeMonteCarlo(sigma,L,d,J,B,T,runs);
        E=energy[Cutoff:runs]
        M=magnetization[Cutoff:runs]

        mean_energy=mean(E)
        Cv=(mean(E.^2)-(mean(E))^2)/T^2
        M2=mean(M.^2)
        Xi=L^d*(mean(M.^2)-(mean(abs.(M)))^2)/T
        
        # Wilson loop stuff is special!
        for k=1:d
            wilson_mean[k]=mean(wilson_loops[Cutoff:runs,k]);
        end


        push!(Temperature,T)
        push!(SpecificHeat,Cv)
        push!(MagnetizationSquared,M2)
        push!(EnergyMean,mean_energy)
        push!(Susceptibility,Xi)
        Wilson_loops_mean=[Wilson_loops_mean;wilson_mean']
        
        next!(p)
    end

return Temperature, SpecificHeat, MagnetizationSquared, EnergyMean, Susceptibility, Wilson_loops_mean

end

