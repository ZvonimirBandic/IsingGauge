
function neighbors(L,d)
    # d-dimensional lattice with length L in each direction
    #
    # returns result array:  neighbor[j,m,1:2]; m runs from 1:d
    #                        provides site index value of neighbors to j in direction "+1" and "-1"
    
    
    N=L^d; # number of sites in the lattice

    # observe that 0:N-1 array would be more convenient than 1:N
    # this brings a lot of -1 and +1 into different array formulas
    
    # mapping of spins from j to (i1,i2,...id) is:
    # j=(id)*L^(d-1)+...+(i_{k-1})*L^(k-2)+...+(i1)+1
	# i_k's are running from 0 to L-1 as it is more convenient for mod function
    # and inverse is:
    # i1=(j-1)mod L
    # ik=round( ((j-1)mod L^k)/L^(k-1))
    # id=round((j-1)/L^(d-1)) # because (any)mod l^d =(any)
       
    # pre-compute neighbors
    # not very optimal, but this computation is done only once
    
    neighbor=zeros(Int32,N,2*d,2) # as there are 2 eighbors in each of the 2 directions
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
        #                        provides site index value of neighbors to j in m-th direction , either in "+1" or "-1"
end
