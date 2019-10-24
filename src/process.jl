## this file runs and processes the main simulations
## it considers the fact that simulations are run parallel.
## it is advised to run initcluster.jl first to setup parallel workers.

# age at infection. 


using Parameters
using Dates

function single(hsize, vc, vs, vt, beta0, beta1, ii, obsize, obtime, hicov, mtime)
    # this function runs 500 simulations with the given parameters.
    # returns a dataframe with averages of 500 simulations
    cd = pmap(1:500) do x
        main(x, hsize, vc, vs, vt, beta0, beta1, ii, obsize, obtime, hicov, mtime)
    end
    
    for i = 1:500
        dts = cd[i]
        insertcols!(dts, 1, :sim => i)    
        insertcols!(dts, 1, :year => 1:mtime)
    end

    d = vcat([cd[i] for i = 1:length(cd)]...)
    ## create yearly average
    dd = d |> @groupby(_.year) |>
        @map({year=key(_), 
              cnt=length(_),
              pop=mean(_.pop),
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

function scenarios()
    # systematically go over the scenarios. 
    
    # vaccine coverage = 0.75, 0.80, 0.85, 0.9, 0.95
    # beta = 0.21 (3.3)
    # delay = short or long. 

    beta = collect(0.1:0.01:0.2)  # beta values 
    hi = collect(0.5:0.05:0.9)    # herd immunity 
    vc = collect(0.7:0.05:0.9)    # vaccine coverage

    ## create the directory for the results.
    dn = "/data/measles/$(Dates.format(Dates.now(), dateformat"mmdd_HHMM"))"
    mkpath("$dn")
    println("created directory: $dn")
       
    ## calculate R0s 
    println("calculating R0s...")
    rzeros= Dict{Float64, Float64}()
    for b in beta
        rd = single(0.0, "fixed", 5, b, 0.0, 0, 0, 0, 0, 2)        
        push!(rzeros, b => rd.inft[1])
    end
    CSV.write("$dn/rzeros.dat", rzeros)

    ## create a progress bar
    total = length(hi)*length(vc)*length(beta)
    p = Progress(total, 0.1, "running scenario")   # minimum update interval: 1 second
    for h in hi, v in vc, b in beta 
        # rdf = single(10000, 0.9, "fixed", 5, 0.21, 0.0, 0, 1, 500, hi, 5)
        # rdd = single(10000, 0.9, "delay", 5, 0.21, 0.0, 0, 1, 500, hi, 5)
        
        # insertcols!(rdf, 3, :R0 => rzeros[b])    
        # insertcols!(rdd, 3, :R0 => rzeros[b])    
        
        # ff, fd = get_save_name(h, v, b)
        # CSV.write("$dn/$ff", rdf)
        # CSV.write("$dn/$fd", rdd)
        sleep(0.01) 
        next!(p; showvalues = [(:hi,h), (:vc,v), (:beta, b)])
    end
end

function get_save_name(h,v,b)
    hh = replace(string(h*100), "0." => "")
    vv = replace(string(v*100), "0." => "")
    bb = replace(string(b), "." => "_")
    
    fixedname = "hi$(hh)_vc$(vv)_beta$(bb)_fixed.dat"
    delayname = "hi$(hh)_vc$(vv)_beta$(bb)_delay.dat"
    return fixedname, delayname    
end




