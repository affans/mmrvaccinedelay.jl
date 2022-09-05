## File is deprecated 
## use scripts/run.jl to run simulations


using Parameters
using Dates
using GLM

function single(hsize, vc, vs, sh, sc, vt, beta0, beta1, obsize, obtime, hicov, mtime)
    # this function runs 500 simulations with the given parameters.
    # returns a dataframe with averages of 500 simulations
    cd = pmap(1:500) do x
        main(x, hsize, vc, vs, sh, sc, vt, beta0, beta1, obsize, obtime, hicov, mtime)
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

function scenarios()
    # systematically go over the scenarios. 
    
    # vaccine coverage = 0.75, 0.80, 0.85, 0.9, 0.95
    # beta = 0.21 (3.3)
    # delay = short or long. 

    #beta = collect(0.1:0.01:0.35)  # beta values
    ros = collect(5:1:30) 
    beta = beta_regression(ros)  ## get the exact beta's from regression analysis
    br = Dict(ros .=> beta)      ## make a dictionary out of it. 

    hi = collect(0.7:0.05:0.9)     # herd immunity 
    vc = collect(0.7:0.01:0.96)    # vaccine coverage

    ## create the directory for the results.
    dn = "/data/measles/$(Dates.format(Dates.now(), dateformat"mmdd_HHMM"))"
    mkpath("$dn")
    println("created directory: $dn")
       
    ## calculate R0s 
    println("calculating R0s...")
    rzeros = calculate_reprodnums(beta)  ## these rzeros should be close to our regressed numberss
    CSV.write("$dn/rzeros.dat", rzeros)

    ## create a progress bar
    total = length(hi)*length(vc)*length(beta)
    p = Progress(total, 0.1, "running scenario")   # minimum update interval: 1 second
    for h in hi, v in vc, r in ros 
        b = br[r]
        rdf =  single(10000, v, "fixed", 1, 1, 1, b, 0.05, 1, [250, 1250], h, 1500)
        rdds = single(10000, v, "delay", 1.7, 14, 1, b, 0.05, 1, [250, 1250], h, 1500)
        rddl = single(10000, v, "delay", 1.7, 29.4, 1, b, 0.05, 1, [250, 1250], h, 1500)

        ff, fds, fdl = get_save_name(h, v, r)
        CSV.write("$dn/$ff", rdf)
        CSV.write("$dn/$fds", rdds)
        CSV.write("$dn/$fdl", rddl)
        
        next!(p; showvalues = [(:hi,h), (:vc,v), (:beta, b)])
    end
end

function scenarios_two()
    ros = [8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25]
    #hi = [0.85] 
    vc = collect(0.45:0.01:0.95)    # vaccine coverage
    beta = beta_regression(ros)  ## get the exact beta's from regression analysis
    br = Dict(ros .=> beta)      ## make a dictionary out of it. 
    
    ## create the directory for the results.
    dn = "/data/measles/$(Dates.format(Dates.now(), dateformat"mmdd_HHMM"))"
    mkpath("$dn")
    println("created directory: $dn")
       
    ## calculate R0s 
    println("calculating R0s...")
    rzeros = calculate_reprodnums(beta)  ## these rzeros should be close to our regressed numberss
    CSV.write("$dn/rzeros.dat", rzeros)

    ## create a progress bar
    total = length(vc)*length(beta)
    p = Progress(total, 0.1, "running scenario")   # minimum update interval: 1 second
    for v in vc, r in ros 
        b = br[r]
        rdf =  single(10000, v, "fixed", 1, 1, 1, b, 0.05, 1, [250, 1250], v, 1500)

        rdd6  = single(10000, v, "delay", 1.7, 14.12, 1, b, 0.05, 1, [250, 1250], v, 1500)
        rdd7  = single(10000, v, "delay", 1.7, 16.47, 1, b, 0.05, 1, [250, 1250], v, 1500)
        rdd8  = single(10000, v, "delay", 1.7, 18.82, 1, b, 0.05, 1, [250, 1250], v, 1500)
        rdd9  = single(10000, v, "delay", 1.7, 21.18, 1, b, 0.05, 1, [250, 1250], v, 1500)
        rdd10 = single(10000, v, "delay", 1.7, 23.53, 1, b, 0.05, 1, [250, 1250], v, 1500)
        rdd11 = single(10000, v, "delay", 1.7, 25.88, 1, b, 0.05, 1, [250, 1250], v, 1500)
        rdd12 = single(10000, v, "delay", 1.7, 28.24, 1, b, 0.05, 1, [250, 1250], v, 1500)
        rdd13 = single(10000, v, "delay", 1.7, 30.59, 1, b, 0.05, 1, [250, 1250], v, 1500)
        rdd14 = single(10000, v, "delay", 1.7, 32.94, 1, b, 0.05, 1, [250, 1250], v, 1500)
        rdd15 = single(10000, v, "delay", 1.7, 35.29, 1, b, 0.05, 1, [250, 1250], v, 1500)
        rdd16 = single(10000, v, "delay", 1.7, 37.65, 1, b, 0.05, 1, [250, 1250], v, 1500)
        rdd17 = single(10000, v, "delay", 1.7, 40.00, 1, b, 0.05, 1, [250, 1250], v, 1500)
        rdd18 = single(10000, v, "delay", 1.7, 42.35, 1, b, 0.05, 1, [250, 1250], v, 1500)
        rdd19 = single(10000, v, "delay", 1.7, 44.71, 1, b, 0.05, 1, [250, 1250], v, 1280)
        rdd20 = single(10000, v, "delay", 1.7, 47.06, 1, b, 0.05, 1, [250, 1250], v, 1280)
        rdd21 = single(10000, v, "delay", 1.7, 49.41, 1, b, 0.05, 1, [250, 1250], v, 1280)
        rdd22 = single(10000, v, "delay", 1.7, 51.76, 1, b, 0.05, 1, [250, 1250], v, 1280)
        rdd23 = single(10000, v, "delay", 1.7, 54.12, 1, b, 0.05, 1, [250, 1250], v, 1280)
        rdd24 = single(10000, v, "delay", 1.7, 56.47, 1, b, 0.05, 1, [250, 1250], v, 1280)


        myfilenames = get_save_name_two(v, v, r)
        #println(myfilenames...)
        CSV.write("$dn/$(myfilenames[1])", rdf)
        CSV.write("$dn/$(myfilenames[2])", rdd6)
        CSV.write("$dn/$(myfilenames[3])", rdd7)
        CSV.write("$dn/$(myfilenames[4])", rdd8)
        CSV.write("$dn/$(myfilenames[5])", rdd9)
        CSV.write("$dn/$(myfilenames[6])", rdd10)
        CSV.write("$dn/$(myfilenames[7])", rdd11)
        CSV.write("$dn/$(myfilenames[8])", rdd12)
        CSV.write("$dn/$(myfilenames[9])", rdd13)
        CSV.write("$dn/$(myfilenames[10])", rdd14)
        CSV.write("$dn/$(myfilenames[11])", rdd15)
        CSV.write("$dn/$(myfilenames[12])", rdd16)
        CSV.write("$dn/$(myfilenames[13])", rdd17)
        CSV.write("$dn/$(myfilenames[14])", rdd18)
        CSV.write("$dn/$(myfilenames[15])", rdd19)    
        CSV.write("$dn/$(myfilenames[16])", rdd20)    
        CSV.write("$dn/$(myfilenames[17])", rdd21)    
        CSV.write("$dn/$(myfilenames[18])", rdd22)    
        CSV.write("$dn/$(myfilenames[19])", rdd23)    
        CSV.write("$dn/$(myfilenames[20])", rdd24)    
    
        
        #println("$dn/hi$(hh)_vc$(vv)_beta$(bb)_delay24.dat")
        next!(p; showvalues = [(:hi,v), (:vc,v), (:beta, b)])
    end
end

function get_save_name_two(h,v,b)
    hh = string(Int(round(h*100)))
    vv = string(Int(round(v*100)))
    bb = string(Int(round(b)))
    
    fixedname = "hi$(hh)_vc$(vv)_beta$(bb)_fixed.dat"
    d6 = "hi$(hh)_vc$(vv)_beta$(bb)_delay6.dat"
    d7 = "hi$(hh)_vc$(vv)_beta$(bb)_delay7.dat"
    d8 = "hi$(hh)_vc$(vv)_beta$(bb)_delay8.dat"
    d9 = "hi$(hh)_vc$(vv)_beta$(bb)_delay9.dat"
    d10 = "hi$(hh)_vc$(vv)_beta$(bb)_delay10.dat"
    d11 = "hi$(hh)_vc$(vv)_beta$(bb)_delay11.dat"
    d12 = "hi$(hh)_vc$(vv)_beta$(bb)_delay12.dat"
    d13 = "hi$(hh)_vc$(vv)_beta$(bb)_delay13.dat"
    d14 = "hi$(hh)_vc$(vv)_beta$(bb)_delay14.dat"
    d15 = "hi$(hh)_vc$(vv)_beta$(bb)_delay15.dat"
    d16 = "hi$(hh)_vc$(vv)_beta$(bb)_delay16.dat"
    d17 = "hi$(hh)_vc$(vv)_beta$(bb)_delay17.dat"
    d18 = "hi$(hh)_vc$(vv)_beta$(bb)_delay18.dat"
    d19 = "hi$(hh)_vc$(vv)_beta$(bb)_delay19.dat"
    d20 = "hi$(hh)_vc$(vv)_beta$(bb)_delay20.dat"
    d21 = "hi$(hh)_vc$(vv)_beta$(bb)_delay21.dat"
    d22 = "hi$(hh)_vc$(vv)_beta$(bb)_delay22.dat"
    d23 = "hi$(hh)_vc$(vv)_beta$(bb)_delay23.dat"
    d24 = "hi$(hh)_vc$(vv)_beta$(bb)_delay24.dat"
    return fixedname, d6, d7, d8, d9, d10, d11, d12, d13, d14, d15, d16, d17, d18, d19, d20, d21, d22, d23, d24
end

function get_save_name(h,v,b)
    hh = string(Int(round(h*100)))
    vv = string(Int(round(v*100)))
    bb = string(Int(round(b)))
    
    fixedname = "hi$(hh)_vc$(vv)_beta$(bb)_fixed.dat"
    delayname = "hi$(hh)_vc$(vv)_beta$(bb)_delays.dat"
    delaynamelong = "hi$(hh)_vc$(vv)_beta$(bb)_delayl.dat"
    return fixedname, delayname, delaynamelong  
end

function calculate_reprodnums(betas)
    ## calculate R0s 
    println("calculating R0s...")
    rzeros = DataFrame(b=betas)
    r0s = map(betas) do b
        rd = single(10000, 1.0, "fixed", 1, 1, 0, b, 0.05, 1, [1], 0, 2) 
        rd.inft[1]
    end
    insertcols!(rzeros, 2, :ro => r0s)
    #CSV.write("$dn/rzeros.dat", rzeros)
    return rzeros
end

function beta_regression(ros)
    ## ros = an array of R0s... this function will calculate and return the corresponding beta value.
    ## however, when running the simulations, the used beta value may not return the same R0 because of stochastcity
    ## use the following beta values to regress. 
    beta = collect(0.1:0.01:0.35)  # beta values
    rs = calculate_reprodnums(beta)
    ols = lm(@formula(ro ~ b), rs)
    c = coef(ols)

    cbs = (ros .- c[1])/c[2]  ## these are the calculated betas
    return cbs 
end

# for i=6:30
#     sc = round(4*i/1.7;digits=2)
#     d = Gamma(1.7, sc)
#     m = round(mean(d);digits=2)
#     println("delay of $i months (dist mean: $m in weeks, $(round((m/4);digits=2)) in months), parameters: $(params(d))")
# end

# function calibrate() 
#     rvals = 5:30 # R values of 5 to 30 
#     beta = beta_regression(rvals)  ## get the exact beta's from regression analysis
#     br = Dict(ros .=> beta)      ## make a dictionary out of it. 
# end

# function beta_regression(ros)
#     ## ros = an array of R0s... this function will calculate and return the corresponding beta value.
#     ## however, when running the simulations, the used beta value may not return the same R0 because of stochastcity
#     ## use the following beta values to regress. 
#     beta = collect(0.1:0.01:0.35)  # beta values
#     rs = calculate_reprodnums(beta)
#     ols = lm(@formula(ro ~ b), rs)
#     c = coef(ols)

#     cbs = (ros .- c[1])/c[2]  ## these are the calculated betas
#     return cbs 
# end
