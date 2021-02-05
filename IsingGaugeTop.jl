include("modules.jl")
include("neighbors.jl")
include("setlattice.jl")
include("compute_energy.jl")
include("hastings.jl")
include("metropolis_ig.jl")
include("wilson.jl")
include("IsingGaugeMonteCarlo.jl")
include("measure.jl")



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
CSV.write(FileName,  DataFrame(Temperature',:auto), header=false)
FileName=string("results/MonteCarlo Filename SpecificHeat N=",N," swees=",runs," ", TimeStamp,".csv")
CSV.write(FileName,  DataFrame(SpecificHeat',:auto), header=false)
FileName=string("results/MonteCarlo Filename Wilson_loops_mean N=",N," sweeps=",runs," ", TimeStamp,".csv")
CSV.write(FileName,  DataFrame(Wilson_loops_mean',:auto), header=false)






