echo "batch run..."

mkdir results # need to have results directory where IsingGaugeTop places CSV results files

# this is an example run

julia IsingGaugeTop.jl 12 3 0 0.5 3 0.05 50000
