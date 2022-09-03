using Revise
using Distributed
using Random
using DelimitedFiles
using Distributions
using Base.Filesystem
using DataFrames
using CSV
using ClusterManagers
using Query
using Statistics
using Dates
using ProgressMeter
using UnPack

using mmrvaccinedelay
const mmr=mmrvaccinedelay

# setup parallel workers
if nprocs() == 1
    if gethostname() == "hpc" 
        println("connecting to hpc")
        using ClusterManagers
        addprocs(SlurmManager(512), N=16, topology=:master_worker, exeflags="--project=.")
    else 
        println("connecting to local")
        addprocs(8, exeflags="--project=.") # for local usage
    end
    # if running Revise, all worker processes will be updated with new code
    @everywhere using mmrvaccinedelay
else 
    println("processors already added")
end

## Simulation Parameters
Base.@kwdef mutable struct SimParameters
    hsize::Int = 10000
    vaccov::Float64 = 0.0       # vaccine coverage
    vacscen::String = "fixed"   # "fixed"=ontime, "delay[x]"=delay with x months 
    vactime::Int64 = 0          # vaccine introduction time
    beta0::Float64 = 0.016       # beta transmission values 
    beta1::Float64 = 0.05        # beta transmission values (controls seasonality)        
    obsize::Int64 = 1            # initial outbreak size
    obtime::Vector{Int64} = [250, 1250]  # times of when outbreak should occur
    hicov::Float64 = 0.0         # population immunity at the start of vaccination
    mtime::Int64 = 1500          # time-span of simulations, in weeks
    rzero::Int64 = 0             # this field is not used... use beta0/beta1 to control rzero
    numofsims::Int64 = 500       # the number of simulations to run. 
end

## declare constants for the scale parameter of the vaccine timing Gamma distribution. Do not change. 
const DELAY_SCALE_PARAM = Dict("fixed"   => 1.25, "delay6"  => 14.12, "delay7"  => 16.47, "delay8"  => 18.82,
                               "delay9"  => 21.18, "delay10" => 23.53, "delay11" => 25.88,"delay12" => 28.24, 
                               "delay13" => 30.59, "delay14" => 32.94, "delay15" => 35.29, "delay16" => 37.65, 
                               "delay17" => 40.00, "delay18" => 42.35, "delay19" => 44.71, "delay20" => 47.06,
                               "delay21" => 49.41, "delay22" => 51.76, "delay23" => 54.12, "delay24" => 56.47)


                               ## This function runs n number of simulations with specified parameters. 
## It returns weekly averages of n simulations.
function run_single(par::SimParameters)
    # extract the variables from par. This is because `SimParameters` is not defined on parallel workers. 
    @unpack hsize, vaccov, vacscen, vactime, beta0, beta1, obsize, obtime, hicov, mtime, rzero, numofsims = par
    delay_scale_param = DELAY_SCALE_PARAM[vacscen]
    (numofsims > 20 && nprocs() < 10) && error("too many simulations, too little processing power")
    cd = pmap(1:numofsims) do x
        main(x, hsize, vaccov, vacscen, 1.7, delay_scale_param, vactime, beta0, beta1, obsize, obtime, hicov, mtime)
    end
    
    for i = 1:numofsims
        dts = cd[i]
        insertcols!(dts, 1, :sim => i)    
        insertcols!(dts, 1, :year => 1:mtime)
    end

    d = vcat([cd[i] for i = 1:length(cd)]...)
    ## create yearly average
    dd = d |> @groupby(_.year) |>
        @map({year=key(_), 
              cnt=length(_),
              avg6=mean(_.avg6),
              susc=mean(_.susc),
              proc=mean(_.proc),
              avg4=mean(_.avg4),
              avg5=mean(_.avg5),
              vinf=mean(_.vinf),
              inft=mean(_.inft), 
              prev=mean(_.prev), 
              reco=mean(_.reco), 
              leftsys=mean(_.leftsys), 
              leftinf=mean(_.leftinf), 
              vacc=mean(_.vacc),
              beta=mean(_.beta),
              meet=mean(_.meet)}) |> DataFrame
    return dd   
end
