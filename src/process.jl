## this file runs and processes the main simulations
## it considers the fact that simulations are run parallel.
## it is advised to run initcluster.jl first to setup parallel workers.

# age at infection. 


using Parameters

function single(vc, vs, vt, beta0, beta1, ii, mtime)
    # this function runs a set of simulations.
    cd = pmap(1:500) do x
        main(x, vc, vs, vt, beta0, beta1, ii, mtime)
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




