function setlattice(L,d,sector)
    
    # this function generates initial desired state of the spin array in the lattice    
    # pick a sector as an input vector sector=(sector_1, sector_2, ..., sector_k,....sector_d) (value +-1)
    # for desired average <W_k> Wilson loop = sector[k]
	# this is how you "store" initial (qu)bits - up to 2^d (qu)bits can be stored
    
    N=L^d # number of lattice sites
    
    #sigma = rand((-1,1),N,d)   #cannot control in advance what values Wilson loops will acquire
    sigma = -ones(N,d) # typical starting state - it does not have electrical charges
	
    # setting lattice default to -1; # implies that initial Wilson loops in all directions are (-1)^L
									 # if sector[k]=1 for all k is chosen
									 
    
    for k=1:d
        if (sector[k] != 1)
            for b=0:(L^(d-k)-1)
                j0=b*L^k
                for s=1:L^(k-1)
                    j1=j0+s  # starting value of j
                    sigma[j1,k]=1 # only first row is set to +1, so all Wilson loops in k direction are "-1"
					# notice that initial value of Wilson loops in k direction becomes (-1)^(L-1)*1
                end
            end
        end
    end
    
    return sigma # returns array of spins
    # sigma[j,m]  j=1:N; d=1:m provides spins on the bonds 1:d around site j
    
end
