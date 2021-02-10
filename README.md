# IsingGauge
Monte Carlo simulator of d-dimensional Ising Gauge hamiltonian

This software is Monte Carlo simulator of Ising Gauge Hamiltonian in arbitrary d dimensions using periodic boundary conditions.
This simulator is especially interesting for studying d=3 case, where classical 3-dimensional Ising Gauge lattice
acts as a topologicaly protected memory, offering storage of 2^3=8 topological (qu)bits (encoded in variable sector).


Top module of the code is called IsingGaugeTop.jl. It includes following modules:

# modules.jl
This modules lists all Julia libraries used in the code


# neighbors.jl
function neighbors(L,d)
L - linear length of array in one direction
d - dimension
N=L^d - total number of spins in the array

We are assuming d-dimensional lattice with length L in each direction. This function returns result array:  neighbor[j,m,1:2]; m runs from 1:d and it provides site index value of neighbors to j in direction "+1" and "-1";

    # mapping of spins from j to (i1,i2,...id) is:
    # j=(id)*L^(d-1)+...+(i_{k-1})*L^(k-2)+...+(i1)+1
	# notice that i_k's are running from 0 to L-1 so that mod functions are easy to calculate
    # and inverse is:
    # i1=(j-1)mod L
    # ik=round( ((j-1)mod L^k)/L^(k-1))
    # id=round((j-1)/L^(d-1)) # because (any)mod l^d =(any)
       
Algorithm for pre-compute of neighbors is not very optimal, but this computation is done only once. This function returns array neighbor[j,m,1:2]; 
where m runs from 1:d. It provides site index value of neighbors to j in direction "+1" and "-1".
    
# setlattice.jl
function setlattice(L,d,sector)

This function acts as a "STORE" and can store 2^d (qu)bits in the array by setting sector[1:d] values to -1 or +1.
This generates initial desired state of the spin array in the lattice    
It allows one to pick a sector as an input vector sector=(sector_1, sector_2, ..., sector_k,....sector_d) (value +-1) for desired average topologically protected <W_k> Wilson loop = sector[k]. It returns sigma which is array of spins:
sigma[j,m]  j=1:N; d=1:m provides spins on the bonds 1:d around site j.

What we expect to see is that for d=2 Wilson loops averages tp stay preserved at T=0 temperature, and for d>2,
we expect Wilson loops averages to stay preserved for temperatures smaller than the phase transition
temperature. For d=3, the phase transition temperature is a continuous phase transition.
   

# compute_energy.jl
function compute_energy(sigma,L,d,J,B,neighbor)
J - Ising gauge strength
B - external field (needed for bias)

This functioncomputes initial energy and magnetization of the lattice, and returns return E, M which are total energy of the lattice and magnetization.

# hastings.jl
function hastings(sigma,L,d,neighbor)
This function performs a Hastings sweep through the lattice. It applies Ising Gauge operator Q,  N times, at random.
Since [H,Q]=0, there should be no energy change (verified). Notice that neighbor comes from caller, to avoid re-computing.
The function returns array sigma, with changes caused by Hastings sweep.
    

# metropolis_ig.jl
function metropolis_ig(E,M,sigma,L,d,J,B,beta,neighbor)
beta=1/T; inverse of temperature and we assume Boltzman constant = 1

Metropolis Monte Carlo single sweep for Ising Gauge lattice and it does flipping bonds degrees of freedom. This function returns sigma,E,M
returns updated energy and magnetization after one monte carlo Metropolis sweep.

# wilson.jl
function wilson(sigma,L,d,neighbor)
    
This function returns averaged Wilson loops across all dimensions k=1:d. All Wilson loops are computed, and then averaged to a single number
The current implementation is debugged, but is not optimal and needs some thoughts whether it is better to explicitely compute neighbors around the loop, or to, compute the coordinates.
    
# IsingGaugeMonteCarlo.jl
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

    # loop for monte carlo sweeps
    
    for s=1:runs
        
        # perform modified Metropolis-Hastings Monte Carlo run
        # Hastings sweep first
        
        sigma,M=hastings(sigma,L,d,neighbor)
        # performs Hastings sweep through the lattice - applies Gauge operator Q, N times
        
        sigma,E,M=metropolis_ig(E,M,sigma,L,d,J,B,beta,neighbor)
        # performs Metropolis sweep through the lattice, and accepts the appropriate flips

        # Wilson loop calculations
        
        W=wilson(sigma,L,d,neighbor) # returns array of averaged Wilson loops across d dimensions
            
    push!(energy,E/N);
    push!(magnetization,M/N);
    wilson_loops=[wilson_loops;W']
       
    end
    

    
    
   return energy, magnetization, sigma, wilson_loops
    
end


# measure.jl
function measure(sector, L,d,J,B,Tmin,Tmax,Tstep,runs

    # average magnetization, average energy and extracting Cv and susceptibility as a function of T
    # Cv=(mean(E^2)-(mean(E))^2)/T^2
    # Xi=L*L*(mean(M.^2)-(mean(abs.(M)))^2)/T 
    # actual susceptibility requires to multiply with number of spins!
    # notice abs(M) - essential for correct behaviour of Xi for scaling, due to SSM happening for real
    # in small lattices
    # sector: indicates which Wilson loop sector is addresssed (+-1 across d dimensions)
    
This module returns measured Temperature, SpecificHeat, MagnetizationSquared, EnergyMean, Susceptibility and averaged Wilson_loops_mean.

# IsingGaugeTop.jl
This is the top Julia module, which is actually being called. It invokes measure function
and has 7  arguments ARGS[1]...ARGS[7] used as:

N=parse(Int,ARGS[1]);   # pick Lmin
d=parse(Int,ARGS[2]);
J=1 # hardcoded
B=parse(Float64,ARGS[3]);
Tmin=parse(Float64,ARGS[4])
Tmax=parse(Float64,ARGS[5])
Tstep=parse(Float64,ARGS[6])
runs=parse(Int,ARGS[7]);

It saves the data in CSV files with a timestamp:

TimeStamp=string(Dates.today()," ", Dates.hour(now())," hr ", Dates.minute(now())," min ", Dates.second(now())," sec")
FileName=string("results/MonteCarlo Filename Temperature N=",N," sweeps=",runs," ", TimeStamp,".csv")
CSV.write(FileName,  DataFrame(Temperature',:auto), header=false)
FileName=string("results/MonteCarlo Filename SpecificHeat N=",N," swees=",runs," ", TimeStamp,".csv")
CSV.write(FileName,  DataFrame(SpecificHeat',:auto), header=false)
FileName=string("results/MonteCarlo Filename Wilson_loops_mean N=",N," sweeps=",runs," ", TimeStamp,".csv")
CSV.write(FileName,  DataFrame(Wilson_loops_mean',:auto), header=false)
