module mmrvaccinedelay

using Parameters      ## with julia 1.1 this is now built in.
using ProgressMeter   ## can now handle parallel with progress_pmap
using DataFrames
using Distributions
using StatsBase
using StaticArrays
using Random

@enum HEALTH UNDEF = 0 SUSC = 1 INF = 2 REC = 3

mutable struct Human ## mutable structs are stored on the heap 
    idx::Int64
    health::HEALTH
    swap::HEALTH   
    infweek::Int64  
    age::Int64          # in weeks 
    ageofdeath::Int64   # in weeks
    group::Int64     
    vaccinated::Bool
    vaccinetime::Int64  # in weeks 
    protection::Float64
    Human() = new()
end

## MAIN SYSTEM PARAMETER
@with_kw mutable struct ModelParameters @deftype Float64
    # general parameters
    sim_time::Int64 = 2500  ## 50 years in 1 week intervals, 50 weeks/year. 
    vaccine_onoff::Bool = false
    vaccination_coverage = 0.8
    vaccination_scenario::String = "fixed"
    herdimmunity_coverage = 0.0
    beta0 = 0.016  ## first results: 0.016
    beta1 = 0.9    ## used for the seasonality amplitude. 
    initial_infected::Int64 = 20
    delay_distribution::Distribution = Gamma(1.7, 14)
end

# global system settings 
# US age distribution for 2018. Downloaded from https://www.census.gov/data/tables/time-series/demo/popest/2010s-national-detail.html 
# the file is on the cluster in /data/datasets/
const agedist =  Categorical(@SVector [0.04826002, 0.061624508, 0.063477388, 0.06418943, 0.066031645, 0.071781172, 0.068554794, 0.066528651, 0.061066338, 0.062600017, 0.06295431, 0.06740441, 0.063575735, 0.054156946, 0.117794637])
const agebraks = @SVector [0:199, 200:449, 450:699, 700:949, 950:1199, 1200:1449, 1450:1699, 1700:1949, 1950:2199, 2200:2449, 2450:2699, 2700:2949, 2950:3199, 3200:3449, 3450:5000]
const humans = Array{Human}(undef, 1000)
const fixed_distribution = Gamma(1.7, 1.25)
const P = ModelParameters()
## code for age bracket generate @SVector [(i):(i+250-1) for i = 200:250:4950]

export HEALTH, humans, agedist, agebraks, P

Base.show(io::IO, ::MIME"text/plain", z::Human) = dump(z)

function main(simnumber=1, hsize=1000, vc=0.8, vs="fixed", shape = 1.7, scale = 29.4, vt=1, b0=0.016, b1=0.9, obsize=1, obtime=[1], hicov=0.0, mtime=2500)  # P is model parameters
    # main entry point of the simulation
    Random.seed!(simnumber)

    # set the global parameters of the model
    P.vaccination_scenario = vs  #shape,scale parameters are not used if scenario == fixed
    P.vaccination_coverage = vc
    P.vaccine_onoff  = false
    P.beta0 = b0
    P.beta1 = b1
    P.sim_time = mtime   
    P.herdimmunity_coverage = hicov
    P.delay_distribution = Gamma(shape, scale)
    # shape 1.7
    # scale=14,  avg delay: 6 months
    # scale=29.4 avg delay: 1 year

    # resize the main humans array
    resize!(humans, hsize)

    ## setup seasonal betas, if P.beta1 = 0 then we have fixed beta
    betas = [P.beta0*(1+P.beta1*sin(2*pi*t/50)) for t = 1:P.sim_time]    
    
    ## data collection variables, might not need all of them. 
    vacc_ctr = zeros(Int64, P.sim_time)   ## counts how many vaccinated at end of week   
    inft_ctr = zeros(Int64, P.sim_time)   ## how many people got infected in a week.    
    prev_ctr = zeros(Int64, P.sim_time)   ## how many people remained infected by end of week (essentiall inft - people who left)
    reco_ctr = zeros(Int64, P.sim_time)   ## how many people recovered by end of week
    left_ctr = zeros(Int64, P.sim_time)   ## how many people left the system by end of week
    left_inf = zeros(Int64, P.sim_time)   ## how many infected people left the system by end of week
    susc_ctr = zeros(Int64, P.sim_time)   ## how many susceptible at each week
    meet_ctr = zeros(Int64, P.sim_time)   ## how many contacts each week
    beta_ctr = zeros(Float64, P.sim_time)
    proc_ctr = zeros(Float64, P.sim_time) 
    avg4_ctr = zeros(Int64, P.sim_time)
    avg5_ctr = zeros(Int64, P.sim_time)
    vinf_ctr = zeros(Int64, P.sim_time)   ## how many infected out of the vaccinated
    avg6_ctr = zeros(Float64, P.sim_time)   ## average age of infection at time t. 

    init_population()    
    init_herdimmunity(P.herdimmunity_coverage)
    
    agm = Vector{Vector{Int64}}(undef, 15) # array of array, holds stratified population

    ## start main time loop    
    for t = 1:P.sim_time         
        ## update agm at every time step
        for grp = 1:15
            agm[grp] = findall(x -> x.group == grp, humans)    
        end 

        ## start vaccineat at specified time with a specified coverage
        if t == vt 
            P.vaccine_onoff = true
            init_vaccination(P.vaccination_coverage)
        end

        ## start a measles outbreak at a specified time with some number
        if t in obtime 
            insert_infected(obsize)
        end

        prev_ctr[t] = length(findall(x -> x.health == INF, humans))  
        avg4_ctr[t] = length(findall(x -> x.health == INF && x.age <= 200, humans)) 
        avg5_ctr[t] = length(findall(x -> x.health == INF && x.age > 200, humans)) 
        vinf_ctr[t] = length(findall(x -> x.health == INF && x.vaccinated == true, humans))
        reco_ctr[t] = length(findall(x -> x.health == REC, humans))
        susc_ctr[t] = length(findall(x -> x.health == SUSC, humans))
        proc_ctr[t] = proc_ctr_update()
        avg6_ctr[t] = avg_age_infection()

        gridsize = length(humans) 
        ## loop for contact tranmission dynamics
        for j in 1:gridsize            
            x = humans[j]
            d, c = contact_dynamic2(x, betas[t], agm) 
            inft_ctr[t] += d
            meet_ctr[t] += c
        end 
        update_swaps()  # move inf -> rec, move susc -> inf if infected.
        
        ## loop for system update 
        for j in 1:gridsize
            x = humans[j]            
            ls, li = age_and_death(x)  ## increase age of each agent. set_vaccine_time if neccessary                        
            vc = vaccinate(x)
            apply_protection(x) # apply their protection level on their new age/vaccination

            left_ctr[t] += ls
            left_inf[t] += li  
            vacc_ctr[t] += vc    
        end        
    end   
    dt = DataFrame(susc=susc_ctr, proc=proc_ctr, inft=inft_ctr, prev=prev_ctr, avg4=avg4_ctr, avg5=avg5_ctr, vinf = vinf_ctr, reco=reco_ctr, leftsys=left_ctr, leftinf=left_inf, vacc=vacc_ctr, beta=betas, meet=meet_ctr, avg6=avg6_ctr)
    return dt
end
export main

function proc_ctr_update()
    # this function calculates the total "amount of protection" 
    # in the population at any given time.
    # used in the main function as a data collection tool 
    protected_humans = findall(x -> x.health == SUSC && x.protection > 0, humans)
    totalprotection = 0.0
    for i in protected_humans
        totalprotection += humans[i].protection
    end
    return totalprotection
end
export proc_ctr_update

function avg_age_infection()
    inf_humans = findall(x -> x.health == INF, humans)
    ages = [humans[i].age for i in inf_humans]
    if length(ages) > 0 
        return mean(ages)
    else 
        return 0 
    end
end
export avg_age_infection

## initialization functions 
function reset_human(x::Human, idx = 0)
    # reset a human back to default values.
    x.idx = idx
    x.health = SUSC
    x.swap = UNDEF 
    x.infweek = 0
    x.group = 1    
    x.age = 0
    x.ageofdeath = 1
    x.vaccinated = false
    x.vaccinetime = 9999
    x.protection = 0
end
export reset_human

function newborn(x::Human)
    # this function resets the human to a default NEW BORN state.
    # we set the age, applying the correct protection, and set a vaccine time is possible.
    reset_human(x, x.idx)
    x.age = 0
    x.group = 1
    x.ageofdeath = calc_ageofdeath(x.age)
    ## give vaccination
    if P.vaccine_onoff 
        if rand() < P.vaccination_coverage           
            set_vaccine_time(x)
        end
    end     
    apply_protection(x) ## should apply maternal immunity
end
export newborn

function init_population()   
    @inbounds for i = 1:length(humans) 
        humans[i] = Human()              ## create an empty human
        reset_human(humans[i], i)        ## reset the human
        apply_agedistribution(humans[i])
        apply_protection(humans[i])
    end
end
export init_population

function init_herdimmunity(cov)
    cov > 1 && error("herd immunity coverage must be between 0 and 1")
    hi_eligble = findall(x -> x.age >= 200 && x.health == SUSC, humans)
    if length(hi_eligble) > 0 
        num = Int(round(cov * length(hi_eligble)))
        h = sample(hi_eligble, num; replace = false)
        @inbounds for i in h 
            humans[i].health = REC
            apply_protection(humans[i])
        end
    end  
end
export init_herdimmunity

function set_vaccine_time(x::Human)
    # we don't set x.vaccinated = true here. age function will take care of it.
    dd = P.delay_distribution
    if P.vaccination_scenario == "fixed"
        x.vaccinetime = 50 + Int(round(rand(fixed_distribution)))
    elseif P.vaccination_scenario == "delay"
        x.vaccinetime = 62 + Int(round(rand(dd)))
    end    
end
export set_vaccine_time

function vaccinate(x::Human)
    ## basic helper function. 
    if x.age == x.vaccinetime 
        x.vaccinated = true
        return true
    end
    return false
end
export vaccinate

function init_vaccination(cov)
    ## everyone between 1 year of age and 4 years of age
    ## will get vaccinated with some efficacy.
    if P.vaccine_onoff
        kids = findall(x -> x.age >= 50 && x.age < 200, humans)
        @inbounds for i in kids
            if rand() < cov     
                x = humans[i]
                set_vaccine_time(x)
                if x.age >= x.vaccinetime 
                    x.vaccinated = true
                end
            end
        end
    end    
end
export init_vaccination

## human functions and helper functions
function apply_agedistribution(x::Human)
    ## for a human, this random assigns an age group, age, and age of death
    x.group = rand(agedist)     
    x.age = rand(agebraks[x.group])
    x.ageofdeath = calc_ageofdeath(x.age)
end
export apply_agedistribution

function insert_infected(num) 
    ## inserts a number of infected people in the population randomly
    l = findall(x -> x.health == SUSC, humans)
    if length(l) > 0 
        h = sample(l, num; replace = false)
        @inbounds for i in h 
            humans[i].health = INF
        end
    end    
end
export insert_infected

function calc_ageofdeath(age)
    ## given a age, what is the chance they will die the year after 
    ## or any of the following years. 
    ## based on USA life tables
    rt = 100 # default age of death
    ageinyears = Int(floor(age/50))    
    for age = ageinyears:100  ## 100% fatal at age 100
        getprobdeath = death_dist[age + 1]  ## since arrays are 1 based
        #println("prob: $getprobdeath")
        if rand() < getprobdeath
            #println(age)
            rt = age
            break
        end
    end
    return rt*50    
end

function get_group(age)
    # input age: in weeks, output group from 1:15 
    r = findfirst(brak -> age in brak, agebraks) 
    if r === nothing # this should never happen
        error("get_group threw an error")
    end
    return r 
end
export get_group

function apply_protection(x)
    # applies a protection level to an individual based on their age/vaccine status
    p = 0.0
    age = x.age
    if x.health == SUSC 
        if age < 26
            p = 0.0
        end    
        if age >= 26 && age < 50
            p = 0.0
        end
        if x.vaccinated
            if P.vaccination_scenario == "fixed"
                p = 0.93
            else 
                p = 0.97
            end 
        end
    end    
    if x.health == REC
        p = 1.0
    end
    x.protection = p
    return p
end
export apply_protection

## transmission dynamics functions
function update_swaps()
    # update swaps, t is current week of the simulation from 1:maxtime
    @inbounds for x in humans
        if x.health == INF
            x.swap = REC
            ## if the person is recovering, the "recovered" population goes up (virtual population)
            ## if the person who is recovering was supposed to die before the simulation ended, 
            ## substract 1 from the virtual population. 
            # popctr[t] = popctr[t] + 1
            # alive_left = x.ageofdeath - x.age # how many weeks remaining till death 
            # if (t+alive_left <= P.sim_time)   # this recovered individual would've died before end of simulation
            #     popctr[t+alive_left] = popctr[t+alive_left] - 1
            # end
        end        
        if x.swap == SUSC
            # person is switching back to susceptible
            # recalculate age of death. their age remains the same.
            x.health = SUSC
            x.swap = UNDEF
            x.ageofdeath = calc_ageofdeath(x.age)
            x.infweek = 0
        end
        if x.swap == INF
            if x.infweek == 0
                x.health = x.swap
                x.swap = UNDEF
            else 
                x.infweek = 0
            end
        end     
        if x.swap == REC
            x.health = x.swap
            x.swap = UNDEF            
        end        
    end  
end
export update_swaps

function age_and_death(x::Human)
    # this function ages the individual by one.
    # also looks at whether the person will be vaccinated. 

    ## data collection from this function
    left_system = false  # true if individual leaves the system by nat death
    left_infect = false  # true if infected individual leaves system by nat death
   
    # store old information and update age
    oldgrp = x.group
    x.age = x.age + 1 
    
    if x.age < x.ageofdeath
        x.group = get_group(x.age)              
    else
        left_system = true
        x.health == INF && (left_infect = true)
        newborn(x)
    end
   
    ## do not change this order as tests depend on it
    return left_system, left_infect
end
export age_and_death

function contact_dynamic2(x, beta, agm)
    ## right away check if x is infected, otherwise no use.
    #c = cts[x.idx, t]
    cnt_infected = 0
    cnt_meet = 0    
    if x.health == INF  ## note: flipping this around to SUSC didnt make a difference to results, only slowed it down
        ig = x.group                # get the infected person's group
        cnt_meet = sum(rand(NB[ig], 7))     # sample the number of contacts     
        # distribute cnt_meet to different groups based on contact matrix. 
        # these are not probabilities, but proportions. be careful. 
        # going from cnt_meet to the gpw array might remove a contact or two due to rounding. 
        gpw = Int.(round.(CM[ig]*cnt_meet)) 
        # let's stratify the human population in it's age groups. 
        # this could be optimized by putting it outside the contact_dynamic2 function and passed in as arguments               
        # enumerate over the 15 groups and randomly select contacts from each group
        for (i, g) in enumerate(gpw)           
            if length(agm[i]) > 0 
                meet = rand(agm[i], g)    # sample 'g' number of people from this group  
                for j in meet             # loop one by one to check disease transmission
                    y = humans[j] 
                    if y.health == SUSC 
                        if rand() < beta*(1 - y.protection)
                            y.swap = INF      # set the swap. 
                            if rand() < 0.5  
                                y.infweek = 0 
                            else 
                                y.infweek = 1
                            end
                            cnt_infected += 1
                        end
                    end           
                end
            end
        end              
        # hp = findall(x -> x.health == SUSC, humans) 
        # meet = rand(hp, cnt_meet)
        # for j in meet             # loop one by one to check disease transmission
        #     y = humans[j] 
        #     if rand() < beta*(1 - y.protection)
        #         y.swap = INF      # set the swap. 
        #         if rand() < 0.5  
        #             y.infweek = 0 
        #         else 
        #             y.infweek = 1
        #         end
        #         cnt_infected += 1
        #     end                     
        # end
    end
    return cnt_infected, cnt_meet
end
export contact_dynamic2

#### DISTRIBUTIONS
function contact_matrix()
    CM = Array{Array{Float64, 1}, 1}(undef, 15)
    CM[1] = [0.2287, 0.1153, 0.0432, 0.0254, 0.0367, 0.0744, 0.1111, 0.0992, 0.0614, 0.0391, 0.0414, 0.0382, 0.032, 0.0234, 0.0305]
    CM[2] = [0.0579, 0.4287, 0.0848, 0.0213, 0.0175, 0.042, 0.0622, 0.0772, 0.0723, 0.0333, 0.0262, 0.0205, 0.0192, 0.0138, 0.0231]
    CM[3] = [0.0148, 0.0616, 0.5082, 0.0885, 0.0219, 0.0212, 0.032, 0.0542, 0.0717, 0.0464, 0.0256, 0.0147, 0.0098, 0.0096, 0.0198]
    CM[4] = [0.01, 0.0202, 0.0758, 0.5002, 0.0804, 0.037, 0.0287, 0.0443, 0.0549, 0.0663, 0.0341, 0.0185, 0.0087, 0.007, 0.0139]
    CM[5] = [0.0191, 0.0155, 0.0203, 0.1107, 0.2624, 0.1342, 0.0932, 0.0668, 0.0604, 0.0754, 0.0616, 0.0334, 0.0157, 0.0092, 0.0221]
    CM[6] = [0.0376, 0.0339, 0.02, 0.0465, 0.126, 0.1847, 0.1238, 0.0895, 0.0742, 0.069, 0.076, 0.0497, 0.0283, 0.0141, 0.0267]
    CM[7] = [0.0554, 0.0623, 0.0444, 0.0362, 0.0536, 0.1071, 0.1592, 0.1372, 0.1001, 0.0677, 0.0586, 0.0426, 0.0328, 0.0175, 0.0253]
    CM[8] = [0.0438, 0.0648, 0.0501, 0.0354, 0.0472, 0.0778, 0.1034, 0.1613, 0.1343, 0.0838, 0.0658, 0.0374, 0.0381, 0.0275, 0.0293]
    CM[9] = [0.0462, 0.0473, 0.0616, 0.0669, 0.0691, 0.0679, 0.0929, 0.1046, 0.1441, 0.1045, 0.0772, 0.0306, 0.0287, 0.0204, 0.0380]
    CM[10] = [0.0237, 0.0254, 0.0409, 0.09, 0.0658, 0.0763, 0.0867, 0.0984, 0.1122, 0.1369, 0.0962, 0.0512, 0.0301, 0.0157, 0.0505]
    CM[11] = [0.0184, 0.0356, 0.0487, 0.0585, 0.0766, 0.098, 0.0804, 0.0811, 0.0922, 0.1018, 0.1109, 0.0818, 0.0439, 0.0183, 0.0538]
    CM[12] = [0.0231, 0.0294, 0.0309, 0.0416, 0.0543, 0.0766, 0.0896, 0.0776, 0.0933, 0.0834, 0.1025, 0.1293, 0.0687, 0.0395, 0.0602]
    CM[13] = [0.0312, 0.0245, 0.0294, 0.0295, 0.0424, 0.0676, 0.0955, 0.0923, 0.0879, 0.0696, 0.0722, 0.1013, 0.1063, 0.0659, 0.0844]
    CM[14] = [0.0212, 0.0357, 0.0251, 0.0234, 0.0346, 0.066, 0.0779, 0.0979, 0.0975, 0.0637, 0.0673, 0.0741, 0.0945, 0.1021, 0.119]
    CM[15] = [0.0202, 0.0276, 0.0508, 0.0539, 0.0315, 0.0513, 0.055, 0.0639, 0.086, 0.089, 0.0677, 0.0594, 0.0755, 0.084, 0.1842]
    return CM
end
const CM = contact_matrix()
export contact_matrix, CM

function negative_binomials()
    ##Age group's mean
    AgeMean = Vector{Float64}(undef, 15)
    AgeSD = Vector{Float64}(undef, 15)

    AgeMean = [10.21, 14.81, 18.22, 17.58, 13.57, 13.57, 14.14, 14.14, 13.83, 13.83, 12.3, 12.3, 9.21, 9.21, 6.89]
    AgeSD = [7.65, 10.09, 12.27, 12.03, 10.6, 10.6, 10.15, 10.15, 10.86, 10.86, 10.23, 10.23, 7.96, 7.96, 5.83]

    # nbinoms = Vector{NegativeBinomial{Float64}}(undef, 15)
    # for i = 1:15
    #     p = 1 - (AgeSD[i]^2-AgeMean[i])/(AgeSD[i]^2)
    #     r = AgeMean[i]^2/(AgeSD[i]^2-AgeMean[i])
    #     nbinoms[i] =  NegativeBinomial(r, p)
    # end
    p = 1 - (AgeSD[1]^2-AgeMean[1])/(AgeSD[1]^2)
    r = AgeMean[1]^2/(AgeSD[1]^2-AgeMean[1])
    nc1 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[2]^2-AgeMean[2])/(AgeSD[2]^2)
    r = AgeMean[2]^2/(AgeSD[2]^2-AgeMean[2])
    nc2 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[3]^2-AgeMean[3])/(AgeSD[3]^2)
    r = AgeMean[3]^2/(AgeSD[3]^2-AgeMean[3])
    nc3 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[4]^2-AgeMean[4])/(AgeSD[4]^2)
    r = AgeMean[4]^2/(AgeSD[4]^2-AgeMean[4])
    nc4 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[5]^2-AgeMean[5])/(AgeSD[5]^2)
    r = AgeMean[5]^2/(AgeSD[5]^2-AgeMean[5])
    nc5 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[6]^2-AgeMean[6])/(AgeSD[6]^2)
    r = AgeMean[6]^2/(AgeSD[6]^2-AgeMean[6])
    nc6 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[7]^2-AgeMean[7])/(AgeSD[7]^2)
    r = AgeMean[7]^2/(AgeSD[7]^2-AgeMean[7])
    nc7 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[8]^2-AgeMean[8])/(AgeSD[8]^2)
    r = AgeMean[8]^2/(AgeSD[8]^2-AgeMean[8])
    nc8 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[9]^2-AgeMean[9])/(AgeSD[9]^2)
    r = AgeMean[9]^2/(AgeSD[9]^2-AgeMean[9])
    nc9 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[10]^2-AgeMean[10])/(AgeSD[10]^2)
    r = AgeMean[10]^2/(AgeSD[10]^2-AgeMean[10])
    nc10 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[11]^2-AgeMean[11])/(AgeSD[11]^2)
    r = AgeMean[11]^2/(AgeSD[11]^2-AgeMean[11])
    nc11 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[12]^2-AgeMean[12])/(AgeSD[12]^2)
    r = AgeMean[12]^2/(AgeSD[12]^2-AgeMean[12])
    nc12 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[13]^2-AgeMean[13])/(AgeSD[13]^2)
    r = AgeMean[13]^2/(AgeSD[13]^2-AgeMean[13])
    nc13 = NegativeBinomial(r, p)


    p = 1 - (AgeSD[14]^2-AgeMean[14])/(AgeSD[14]^2)
    r = AgeMean[14]^2/(AgeSD[14]^2-AgeMean[14])
    nc14 = NegativeBinomial(r, p)

    p = 1 - (AgeSD[15]^2-AgeMean[15])/(AgeSD[15]^2)
    r = AgeMean[15]^2/(AgeSD[15]^2-AgeMean[15])
    nc15 = NegativeBinomial(r, p)
    return nc1,nc2,nc3,nc4,nc5,nc6,nc7,nc8,nc9,nc10,nc11,nc12,nc13,nc14,nc15
end
const NB = negative_binomials()
export negative_binomials, NB

## see excel file in measles folder for death distribution, this is per year
const death_dist = [0.006799, 0.000483, 0.000297, 0.000224, 0.000188, 0.000171, 0.000161, 0.000151, 0.000136, 
        0.000119, 0.000106, 0.000112, 0.000149, 0.000227, 0.000337, 0.00046, 0.000579, 0.000684, 
        0.000763, 0.000819, 0.000873, 0.000926, 0.00096, 0.000972, 0.000969, 0.00096, 0.000954, 
        0.000952, 0.000958, 0.000973, 0.000994, 0.001023, 0.001063, 0.001119, 0.001192, 0.001275, 
        0.001373, 0.001493, 0.001634, 0.001788, 0.001945, 0.002107, 0.002287, 0.002494, 0.002727,
        0.002982, 0.003246, 0.00352, 0.003799, 0.004088, 0.004404, 0.00475, 0.005113, 0.005488, 
        0.005879, 0.006295, 0.006754, 0.00728, 0.007903, 0.008633, 0.009493, 0.010449, 0.011447, 
        0.012428, 0.013408, 0.014473, 0.015703, 0.017081, 0.018623, 0.020322, 0.022104, 0.024023, 
        0.026216, 0.028745, 0.031561, 0.034427, 0.037379, 0.040756, 0.044764, 0.049395, 0.054471, 
        0.059772, 0.065438, 0.071598, 0.078516, 0.085898, 0.093895, 0.102542, 0.111875, 0.121928, 
        0.132733, 0.144318, 0.156707, 0.169922, 0.183975, 0.198875, 0.21462, 0.231201, 0.2486, 0.266786, 1.0] 
export death_dist

end # module
