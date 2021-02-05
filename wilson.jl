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
