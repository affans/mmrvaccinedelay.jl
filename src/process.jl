## this file runs and processes the main simulations
## it considers the fact that simulations are run parallel.
## it is advised to run initcluster.jl first to setup parallel workers.

# age at infection. 


using Parameters

function single(beta0, beta1, mtime, ii, hic, vc)
    # this function runs a set of simulations.
    cd = pmap(1:500) do x
        main(x, hic, vc, beta0, beta1, ii, mtime)
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
              avg_susc=mean(_.susc),
              avg_inft=mean(_.inft), 
              avg_prev=mean(_.prev), 
              avg_reco=mean(_.reco), 
              avg_leftsys=mean(_.leftsys), 
              avg_leftinf=mean(_.leftinf), 
              avg_vacc=mean(_.vacc)}) |> DataFrame
    return dd
end