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

