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