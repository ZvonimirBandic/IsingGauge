using Plots, Statistics
using ProgressMeter
using LinearAlgebra,    LsqFit;
using SparseArrays,  KrylovKit
using Profile
using DataFrames, CSV, Dates

function neighbors(L,d)
    # d-dimensional lattice with length L in each direction
    #
    # returns result array:  neighbor[j,m,1:2]; m runs from 1:d
    #                        provides site index value of neighbors to j in direction "+1" and "-1"
    
    
    N=L^d; # number of sites in the lattice

    # observe that 0:N-1 array would be more convenient than 1:N
    # this brings a lot of -1 and +1 into different array formulas
    
    # mapping of spins from j to (i1,i2,...id) is:
    # j=(id-1)*L^(d-1)+...+(i_{k-1}-1)*L^(k-2)+...+(i1-1)+1
    # and inverse is:
    # i1=(j-1)mod L
    # ik=round( ((j-1)mod L^k)/L^(k-1))
    # id=round((j-1)/L^(d-1)) # because (any)mod l^d =(any)
       
    # pre-compute neighbors
    # not very optimal, but this computation is done only once
    
    neighbor=zeros(Int32,N,2*d,2) # as there are 2 neighbors in each of the 2 directions
    i=zeros(Int32,d)
    ir=zeros(Int32,d,d) # ir[m,k] is k-th coordinate of the right neighbor in the m-th direction
    il=zeros(Int32,d,d) # il[m,k] is k-th coordinate of the left neighbor in the m-th direction
    
    for j=1:N
        # go through all spins
        
        for k=1:d # for all coordinates k=1 is x, k=2 is y, k=3 is z etc.
            i[k]=floor(mod(j-1,L^k)/L^(k-1)) # coordinates for spin (j-1)
        end
        
        # each spin has 2 neighbors in one of the d directions - total of 2d neighbors
        # for each we need d coordinates
        
        for m=1:d # for each of the directions
            ir[m,1:d]=i[1:d]; # copy coordinates of j into ir
            il[m,1:d]=i[1:d]; # just copy coordinates of j into il
        end
        
        for m=1:d # for each of the directions
            ir[m,m]=mod(ir[m,m]+1,L) # "+1" neighbor in the m-th direction has m-th coordinate impacted!
            il[m,m]=mod(il[m,m]-1,L) # "-1" neighbor in the m-th direction has m-th coordinate impacted
            
        end
        
        
        for m=1:d # directions
            for k=1:d # coordinates
                neighbor[j,m,1]=neighbor[j,m,1]+ir[m,k]*L^(k-1) # right neigbor into m-th direction
                neighbor[j,m,2]=neighbor[j,m,2]+il[m,k]*L^(k-1)
            end
            neighbor[j,m,1]=neighbor[j,m,1]+1;
            neighbor[j,m,2]=neighbor[j,m,2]+1;
        end
        

    end
    
    return neighbor
        # returns result array:  neighbor[j,m,1:2]; m runs from 1:d
        #                        provides site index value of neighbors to j in direction "+1" and "-1"
end

function setlattice(L,d,sector)
    
    # this function generates initial desired state of the spin array in the lattice    
    # pick a sector as an input vector sector=(sector_1, sector_2, ..., sector_k,....sector_d) (value +-1)
    # for desired average <W_k> Wilson loop = sector[k]
    
    N=L^d # number of lattice sites
    
    #sigma = rand((-1,1),N,d)   #This is your inital random state    
    sigma = -ones(N,d) # test only, but worth more investigation - does not have charges at the start!
    # setting lattice default to -1; be very careful with lattices where L=odd number
    # this effect alone is worth further investigation
    
    for k=1:d
        if (sector[k] != 1)
            for b=0:(L^(d-k)-1)
                j0=b*L^k
                for s=1:L^(k-1)
                    j1=j0+s  # starting value of j
                    sigma[j1,k]=1 # only first row is set to -1
                end
            end
        end
    end
    
    return sigma # returns array of spins
    # sigma[j,m]  j=1:N; d=1:m provides spins on the bonds 1:d around site j
    
end


function compute_energy(sigma,L,d,J,B,neighbor)
    
    # compute initial energy and magnetization of the lattice
    # runs through the whole lattice so is O(N)
    # neighbor comes from other functions, to avoid computing multiple times
    # I wander if I should make neighbor global
    
    E=0; M=0;
    N=L^d
    

    
    for j=1:N
        # there should be (d choose 2) plaquettes associated with each spin
        for m=1:d-1
            for n=(m+1):d
                # compute energy of the corresponding plaquette
                # find indices of plaquette spins first
                j1=j; # this is the first corner of the plaquette
                j2=neighbor[j1,m,1]; # neighbor to j1, to the "+1"into m-th direction
                j3=neighbor[j2,n,1]; # neighbor to j2, to the "+1" into the n-th direction
                j4=neighbor[j1,n,1]; # neighbor to the j, to the "+1" into the n-th direction
                j4_prime=neighbor[j3,m,2]; # this needs to be the same as j4
                
                #println("for spin j=",j," the (",m,",",n,")-th plaquette is ", j1,",",j2,",",j3,",",j4)
                
                if (j4 != j4_prime)
                    println("FAILURE")
                end
                
                E=E-J*sigma[j1,m]*sigma[j2,n]*sigma[j4,m]*sigma[j1,n]
            
            end
        end
    end
    
    M=sum(sigma);
    E=E-B*sum(sigma); # this is the term that will favor all spins going "+1"
    return E, M
    # returns total energy of the lattice and magnetization
end


function hastings(sigma,L,d,neighbor)
    # performs a Hastings sweep through the lattice
    # it applies Ising Gauge operator Q,  N times, at random
    # since [H,Q]=0, there should be no energy change (verified)
    # neighbor comes from caller, to avoid re-computing
    
    N=L^d
 
    
    for g=1:N # flip local gauge N times randomly
        j=rand(1:N)
        for m=1:d # looking into d directions
            # local Gauge invariance operator
            sigma[j,m]=-sigma[j,m]
            sigma[neighbor[j,m,2],m]=-sigma[neighbor[j,m,2],m]        
        end
    end
    magnetization=sum(sigma)
    
    return sigma,magnetization
    # returns usual sigma[j,m] spin array and magnetization
end


function metropolis_ig(E,M,sigma,L,d,J,B,beta,neighbor)
    # Metropolis Monte Carlo single sweep for Ising Gauge lattice
    # focuses on flipping bonds degrees of freedom
    
    N=L^d; # total number of spins

    for l=1:N*d # one sweep of M.C. is size of N*d - number of actual spins
            dE=0;dM=0;
            
            # flip a spin
            j=rand(1:N)
            m=rand(1:d)
            
            #println("______________________________________________________")
            #println("test flipping (",j,",",m," )")
            
            # evaluate change in energy and magnetization
            # impacts all plaquettes containing spin sigma(j,m), which I think is (d-1)*2 plaquettes
            
            for n=1:d # explore all directions except one defined by m
                if (n!=m)
                    # now we have j, m and n
                    # this is the plaquette in the "+1" direction along n-th direction
                        j1=j; # this is the first corner of the plaquette
                        j2=neighbor[j1,m,1]; # neighbor to j1, to the "+1"into m-th direction
                        j3=neighbor[j2,n,1]; # neighbor to j2, to the "+1" into the n-th direction
                        j4=neighbor[j1,n,1]; # neighbor to the j, to the "+1" into the n-th direction
                        dE=dE+2*J*sigma[j1,m]*sigma[j2,n]*sigma[j4,m]*sigma[j1,n]
            #println("Neigboring plaquettes are ",j1,",",j2,",",j3,",",j4)
            #println("Spins on this plaquette are: ",sigma[j1,m],", ",sigma[j2,n],", ",sigma[j4,m]," and ",sigma[j1,n])
                    # this is the plaquette in the "-1" direction along the n-th direction
                        j1=j; # this is the first corner of the plaquette
                        j2=neighbor[j1,m,1]; # neighbor to j1, to the "+1"into m-th direction
                        j3=neighbor[j2,n,2]; # neighbor to j2, to the "-1" into the n-th direction
                        j4=neighbor[j1,n,2]; # neighbor to the j, to the "-1" into the n-th direction
                        dE=dE+2*J*sigma[j1,m]*sigma[j3,n]*sigma[j4,m]*sigma[j4,n]
            #println("                      and ",j1,",",j2,",",j3,",",j4)
            #println("and on this plaquette spins are: ",sigma[j1,m],", ",sigma[j3,n],", ",sigma[j4,m]," and ",sigma[j4,n])
                end
            end
            
            dM=-2*sigma[j,m]
            dE=dE-B*dM; # I think I got the sign right, as E=-J*Sum(plaquette terms)-B*Sum(sigma)
            #println("dE= ", dE)
            #println("dM= ", dM)
        
                    # make a decision whether to accept it
            
            if (dE < 0) 
                sigma[j,m]=-sigma[j,m]
                E=E+dE
                M=M+dM
                #println("Flip at (",j,",",m," ) ACCEPTED")
                #println("New energy E= ", E)

            elseif (dE>0)
                r=rand(Float64) # generate random number between 0 and 1
                if (exp(-beta*dE)>r) # MOnte Carlo approach - accept within exp(-dE/(kB*T)
                    sigma[j,m]=-sigma[j,m]
                    E=E+dE
                    M=M+dM
                    #println("r=",r," was small enough, so flip at (",j,",",m," ) ACCEPTED")
                    #println("New energy E= ", E)
                end    
            end
        
    end

    return sigma,E,M
    # returns change in energy and magnetization after one monte carlo Metropolis sweep

end

function wilson(sigma,L,d,neighbor)
    
    # this function returns averaged Wilson loops across all dimensions k=1:d
    # all Wilson loops are computed, and then averaged to a single number
    # the current implementation is debugged, but is not optimal
    # needs some thoughts whether it is better to explicitely compute
    # neighbors around the loop, or to, compute the coordinates
    
    W=ones(Float32,d)    # array that contains final averaged Wilson loops
        
    for k=1:d # we will compute averages in k directions for Wilson loops
        # each direction will have computed L^(k-1) different loops    
        # compute j0 that belongs to the plane defined by i[k]=0
        # then move in steps of j0, j0+L^(k-1), j0+2*L^(k-1) - i.e. neighbors in kth direction
        # I think using neighbor function would likely be faster
            
        wilson_sum=0;
            
        for b=0:(L^(d-k)-1)
            j0=b*L^k
                for s=1:L^(k-1)
                #   println("b= ",b,", s=", s," ")    
                    j1=j0+s  # starting value of j    
                 #   println("j1= ",j1)
                    wilson=1;
                    step=L^(k-1);
                
                    for m=0:(L-1)  # now moving along actual Wilson loop
                        wilson=wilson*sigma[j1,k]                       
                       # println("sigma(",j0+m*L^(k-1)," , ",k,")=   ",sigma[j0+m*L^(k-1),k])
                       # println("sigma(",j1,",",k,")= ",sigma[j1,k])
                        j1=j1+step;  
                    end
                  #  println("wilson loop for j1 =",j1-step , " is =  ", wilson)
                    wilson_sum=wilson_sum+wilson
                end
        end
        W[k]=wilson_sum/(L^(d-1)) # compute average
    end
    return W
    # W contains array of averaged Wilson loops
    
end

function IsingGaugeMonteCarlo(sigma,L,d,J,B,T,runs)
    
    
    # classical (Ising) Gauge Hamiltonian Monte Carlo simulation
    # this is top function, that uses following functions:
    #
    #     neighbors, compute_energy, hastings, metropolis_ig
    #     
    #
    # desired sector - vector of initial state of Wilson loops in T~0 state
    # sigma[j,m] contains initial spins state of the lattice
    # j runds from 1:N, and m runs from 1:d and sigma[j,m] belongs to m bonds connected to site j in "+1" direction
    # L = Linear length in one e_i direction (i=1,2,....d)
    # d - number of dimensions
    # J - energy; H=-J*sigma over plaquettes
    # T - temperature
    
    N=L^d # number of spins
     # we have total of dN spins for d-dimensional Ising gauge H
    beta=1/T

    
    # observe that 0:N-1 array would be more convenient than 1:N
    # this brings a lot of -1 and +1 into different array formulas
    
    # mapping of spins from j to (i1,i2,...id) is:
    # j=(id-1)*L^(d-1)+...+(i_{k-1}-1)*L^(k-2)+...+(i1-1)+1
    # and inverse is:
    # i1=(j-1)mod L
    # ik=round( ((j-1)mod L^k)/L^(k-1))
    # id=round((j-1)/L^(d-1)) # because (any)mod l^d =(any)
    
    # pre-compute neighbors
    # not very optimal, but this computation is done only once
    
    neighbor=neighbors(L,d) 
    # call to function that computes data array neighbor
    # returns result array:  neighbor[j,m,1:2]; m runs from 1:d
    #                        provides site index value of neighbors to j in direction "+1" and "-1"
       

    energy        = Float32[] # array that will store energy computed in every sweep
    magnetization = Float32[] # array that will store magnetization computed in every sweep
    wilson_loops=Array{Int32,2}(undef,0,d) # array of wilson loops in directions 1,2,...d at each temperature step
    
    # compute initial energy and magnetization of the lattice
    # this is computation based on compute_energy function that runs through whole lattice and is slow
    
    E,M=compute_energy(sigma,L,d,J,B,neighbor)
    
    #println("Initial array energy is E= ", E)
    #println("Initial sigma is", sigma)
    
    
    # loop for monte carlo sweeps
    
    for s=1:runs
        #println("this is run= ",s, "out of ",runs)
        #println("E= ", E)
        
        # perform modified Metropolis-Hastings Monte Carlo run
        # Hastings sweep first
        
        sigma,M=hastings(sigma,L,d,neighbor)
        # performs Hastings sweep through the lattice - applies Gauge operator Q, N times
        
        sigma,E,M=metropolis_ig(E,M,sigma,L,d,J,B,beta,neighbor)
        # performs Metropolis sweep through the lattice, and accepts the appropriate flips

         #   println("spins direction -> 1 ",sigma[1:16,1]);
         #   println("spins direction |v 2 ",sigma[1:16,2])
        
        # Wilson loop calculations
        
        W=wilson(sigma,L,d,neighbor) # returns array of averaged Wilson loops across d dimensions
            
    push!(energy,E/N);
    push!(magnetization,M/N);
    wilson_loops=[wilson_loops;W']
       
    end
    

    
    
   return energy, magnetization, sigma, wilson_loops
    
end

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

sector=[1;-1;1]; #hardcoded
N=parse(Int,ARGS[1]);   # pick Lmin
d=parse(Int,ARGS[2]);
J=1 # hardcoded
B=parse(Float64,ARGS[3]);
Tmin=parse(Float64,ARGS[4])
Tmax=parse(Float64,ARGS[5])
Tstep=parse(Float64,ARGS[6])
runs=parse(Int,ARGS[7]);


Temperature,SpecificHeat,MagnetizationSquared,EnergyMean,Susceptibility,Wilson_loops_mean=measure(sector,N,d,J,B,Tmin,Tmax,Tstep,runs)


# saving the data
TimeStamp=string(Dates.today()," ", Dates.hour(now())," hr ", Dates.minute(now())," min ", Dates.second(now())," sec")
FileName=string("results/MonteCarlo Filename Temperature N=",N," sweeps=",runs," ", TimeStamp,".csv")
CSV.write(FileName,  DataFrame(Temperature',:auto), writeheader=false)
FileName=string("results/MonteCarlo Filename SpecificHeat N=",N," swees=",runs," ", TimeStamp,".csv")
CSV.write(FileName,  DataFrame(SpecificHeat',:auto), writeheader=false)
FileName=string("results/MonteCarlo Filename Wilson_loops_mean N=",N," sweeps=",runs," ", TimeStamp,".csv")
CSV.write(FileName,  DataFrame(Wilson_loops_mean',:auto), writeheader=false)






