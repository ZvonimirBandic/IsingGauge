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

