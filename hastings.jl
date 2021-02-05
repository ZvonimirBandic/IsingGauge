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
