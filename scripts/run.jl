using Revise
using Distributed
using ClusterManagers
using Base.Filesystem
using Query
using Statistics
using UnPack
using Gnuplot
using GLM
using DataFrames
using ProgressMeter
using Dates
using DrWatson
using CSV

using mmrvaccinedelay
#const mmr=mmrvaccinedelay

ENV["JULIA_WORKER_TIMEOUT"] = 120  ## julia 1.8 is taking longer to connect to workers

# setup parallel workers
if nprocs() == 1
    if gethostname() == "hpc" 
        println("connecting to hpc, 2 nodes onl, edited file")
        addprocs(SlurmManager(320), N=10, verbose="", topology=:master_worker, exeflags="--project=$(Base.active_project())")
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
    inftime_prob::Float64 = 0.5  # probability of infectiousness in week 1 or week 2 (0.0 = week 1, 1.0 = week 2)
    mtime::Int64 = 1500          # time-span of simulations, in weeks
    rzero::Int64 = 0             # this field is not used... use beta0/beta1 to control rzero
    numofsims::Int64 = 500       # the number of simulations to run. 
end

function run_single(par::SimParameters)
    ## This function runs n number of simulations with specified parameters. 
    ## It returns weekly averages of n simulations.
    ## extract the variables from par. This is because `SimParameters` is not defined on parallel workers. 
    @unpack hsize, vaccov, vacscen, vactime, beta0, beta1, obsize, obtime, hicov, inftime_prob, mtime, rzero, numofsims = par
    
    ## declare constants for the scale parameter of the vaccine timing Gamma distribution. Do not change. 
    DELAY_SCALE_PARAM = Dict("fixed"   => 1.25, "delay6"  => 14.12, "delay7"  => 16.47, "delay8"  => 18.82,
    "delay9"  => 21.18, "delay10" => 23.53, "delay11" => 25.88,"delay12" => 28.24, 
    "delay13" => 30.59, "delay14" => 32.94, "delay15" => 35.29, "delay16" => 37.65, 
    "delay17" => 40.00, "delay18" => 42.35, "delay19" => 44.71, "delay20" => 47.06,
    "delay21" => 49.41, "delay22" => 51.76, "delay23" => 54.12, "delay24" => 56.47)
    delay_scale_param = DELAY_SCALE_PARAM[vacscen]

    (numofsims > 20 && nprocs() < 10) && error("too many simulations, too little processing power")
    #println("starting single run pmap")
    cd = pmap(1:numofsims) do x
        main(x, hsize, vaccov, vacscen, 1.7, delay_scale_param, vactime, beta0, beta1, obsize, obtime, hicov, inftime_prob, mtime)
    end
    
    for i = 1:numofsims
        dts = cd[i]
        insertcols!(dts, 1, :sim => i)    
        insertcols!(dts, 1, :week => 1:mtime)
    end

    d = vcat([cd[i] for i = 1:length(cd)]...)
    ## create weekly average of numofsims simulations
    dd = d |> @groupby(_.week) |>
        @map({week=key(_), 
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

function run_scenarios()
    ## This function sets up the scenarios/parameters 
    ## and calls run_single() for each set. 
    ## Depending on the number of processors and number of scenarios,
    ## this function can take a long time to run. 

    ## create the directory for the results.
    dn = "/data/measles/$(Dates.format(Dates.now(), dateformat"mmdd_HHMM"))"
    mkpath("$dn"); println("saving simulation results to directory: $dn...")
      
    
    ## setup default simulation parameter values
    par = SimParameters() # keep most things default 
    
    # systematically go over the scenarios.     
    br = Dict([12, 14, 16, 18]  .=> get_beta.([12, 14, 16, 18])) # Beta values from the regression analysis
    vc = collect(0.60:0.01:0.95) # vaccine coverage from 60% to 95%.

    total = length(vc)*length(br)

    ## print debug information
    println("running on $(nprocs()) processors...")
    println("default parameter values..."); dump(par)
    println("starting simulations, total of $total scenarios...")

    ## create a progress bar for the simulations
    p = Progress(total, dt=0.1, barglyphs=BarGlyphs("[=> ]"))  
        
    for v in vc, (r, b) in br 
        par.vactime = 1 # all vaccination scenarios begin on day 1.
        par.vaccov = v  # vaccine coverage
        par.beta0 = b
        par.hicov = v   # initial herd immunity  
        par.rzero = r   # r-zero value, not really used
        
        # run on-time "fixed" scenario
        rdf = run_single(par)        
        fn = savename(par, "dat", accesses=[:vaccov, :rzero, :hicov, :vacscen])
        CSV.write("$dn/$fn", rdf)
        
        ## run delay scenarios
        for dsc in ("delay$x" for x = 6:24)
            par.vacscen = dsc
            rdf = run_single(par)
            fn = savename(par, "dat", accesses=[:vaccov, :rzero, :hicov, :vacscen])
            CSV.write("$dn/$fn", rdf)
        end
        next!(p; showvalues = [(:vaccov, par.vaccov), (:beta, par.beta0)])
    end
end


function rvalue_linreg(; plot_figures = false)
    ## input betas: a vector of beta values 
    println("Estimating model R0...")
    

    beta_values = collect(0.1:0.01:0.35)  # beta values
    
    # set up parameters 
    sp = SimParameters() 
    sp.obsize = 1         # initial outbreak size
    sp.obtime = [1]       # times of when outbreak should occur
    sp.mtime = 1    # time-span of simulations, in weeks (everything should happen in 1 week, but run for 2 weeks to get all the other data collection vectors)
    sp.inftime_prob = 1.0  # will force all infected individuals to become infectious in one week. 
    rvals = map(beta_values) do b
        sp.beta0 = b
        rd = run_single(sp)
        rd.inft[1]
    end

    # apply linear regression
    df = DataFrame(b = beta_values, r = rvals)
    ols = lm(@formula(r ~ b), df)
    Yp = predict(ols); #load predictions into Yp
    
    c = coef(ols)
    println("the coefficients of the regression are: $c")
    
    #@gp "reset"
    #@gp :- "set xlabel 'beta value'"
    #@gp :- "set ylabel 'R'"
    #@gp :- beta_values r0s "with lines title 'data' lw 2 lc 'blue'"
    #@gp :- beta_values Yp "with lines title 'line of best fit' lw 3 lc 'black'"
    #display(@gp)
    #Gnuplot.save(term="svg enhanced font 'Arial,10' size 600,250", output="calibration_reg.svg")

    #c = coef(ols)
    #cbs = (ros .- c[1])/c[2]  ## these are the calculated betas
    return DataFrame( (; rvals, beta_values) )
end

function get_beta(r_val) 
    # function returns the calibrated beta value given a reproduction number 
    # this is done through solving the linear regression R = c0 + c1 * beta 
    # where c0 and c1 are determined through the regression in rvalue_linreg
    c0 = -0.19415384615384193
    c1 = 92.41675213675214
    cbs = (r_val .- c0)/c1  ## these are the calculated betas
end