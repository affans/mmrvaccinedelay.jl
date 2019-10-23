using mmrvaccinedelay
using Test

const mmr = mmrvaccinedelay
const gridsize = length(mmr.humans)
const gridtime = 50*50
const vax_coverage = 0.8

# to do:
# test whether group changes when individual is at the end of their age bracket
@testset "INIT" begin
    mmr.init_population()
    @test length(findall(x -> x == undef, mmr.humans)) == 0 ## check for proper initialization
    for i = 1:gridsize
        x = mmr.humans[1]
        @test x.idx == 1
        @test x.vaccinetime == 9999
        @test x.vaccinated == false
        @test x.swap == mmr.UNDEF
        @test x.health == mmr.SUSC
        @test x.protection == 0
        @test x.age <= x.ageofdeath
        @test x.protection == 0 ## no one is vaccinated or protected. 
    end
    
    ## test if the insert_infected function really works.
    mmr.init_population()
    mmr.insert_infected(1)
    a = findall(x -> x.health == mmr.INF, mmr.humans)
    @test length(a) == 1

    ## test if the insert_infected dosn't sample the sample person twice.
    mmr.init_population()
    mmr.insert_infected(100)
    a = findall(x -> x.health == mmr.INF, mmr.humans)
    @test length(a) == 100

    ## test init_herdimmunity function
    mmr.init_population()
    mmr.init_herdimmunity(0)    
    a3 = findall(x -> x.health == mmr.REC, mmr.humans)
    @test length(a3) == 0 

    hicov = 0.8
    mmr.init_population()
    hi_eligble = findall(x -> x.age >= 200 && x.health == mmr.SUSC, mmr.humans)
    mmr.init_herdimmunity(hicov)    
    a3 = findall(x -> x.health == mmr.REC, mmr.humans)
    @test isapprox(length(a3)/length(hi_eligble), hicov; atol=0.1)

    # test if susceptibles are still with 0 protection
    for x in mmr.humans 
        @test x.health == mmr.SUSC || x.health == mmr.REC ## should not get an infected
        @test x.vaccinated == false 
        @test x.vaccinetime == 9999
        if x.health == mmr.SUSC 
            @test x.protection == 0
        elseif x.health == mmr.REC
            @test x.protection == 1
        end
    end
end

@testset "AGE" begin
    mmr.init_population()    
    ## check if age_and_vaccinate function (age only)
    ## checks if age is incremented properly
    ## and verifies dead people are replaced properly. 
    a1 = [mmr.humans[i].age for i = 1:gridsize]
    leave = zeros(Bool, gridsize) #contains whether they have left the system
    for i = 1:gridsize
        ls, li = mmr.age_and_death(mmr.humans[i])
        leave[i] = ls 
    end
    ## check if the age_and_death worked properly. 
    a2 = [mmr.humans[i].age for i = 1:gridsize]
    for i = 1:gridsize
        @test (a2[i] == a1[i] + 1) || a2[i] == 0
    end
    # test whether the number of newborns = number of people left
    l1 = length(findall(x -> x == 0, a2))
    l2 = length(findall(x -> x == true, leave))
    @test l1 == l2

    ## this is not a test. basically runs the model for gridtime and 
    ## stores who left the system. this needs to be plotted or something
    mmr.init_population()
    wd = zeros(Int64, gridtime)
    for t = 1:gridtime
        for x in mmr.humans 
            r = mmr.age_and_death(x)
            if r[2] == true 
                wd[t] += 1
            end
        end
    end
    #UnicodePlots.scatterplot(wd)
end

@testset "SWAPS" begin
    ## This tests whether the person's swaps are set properly.
    ## secondary tests makes sure if person isn't infected, they can't transmit the disease even with beta = 1
    mmr.init_population()
    x = mmr.humans[1]            # find a human
    x.swap = mmr.INF             # set SWAP to infected
    @test x.infweek == 0
    mmr.update_swaps()           # update swap
    @test x.swap == mmr.UNDEF    # tests whether swap happened
    @test x.health == mmr.INF
 
    # go through one week and see what happens after 
    agm = Vector{Vector{Int64}}(undef, 15) # array of array, holds stratified population
    for grp = 1:15
        agm[grp] = findall(x -> x.health == mmr.SUSC && x.group == grp, mmr.humans)    
    end 
    for i in 1:gridsize
        d, c = mmr.contact_dynamic2(mmr.humans[i], 1.0, agm)
        # we only made the first human infected, so the number of secondary infections should be zero for everyone else        
    end 
    @test x.health == mmr.INF    # human[1] is still infected 
    mmr.update_swaps()           # this should move their health to REC
    @test x.health == mmr.REC    # test whether it recovers back to REC (it will fail if INC -> SUSC depending on vaccine/prevaccine era)
    @test x.swap == mmr.UNDEF
    # test the population to see only one recovered
    @test 1 == length(findall(x -> x.health == mmr.REC, mmr.humans))
    @test length(findall(x -> x.health == mmr.INF, mmr.humans)) >= 0 
end

@testset "VAX" begin
    ## test vaccination coverage and protection level
    mmr.init_population()
    mmr.init_vaccination(vax_coverage)
    v = 0
    for x in mmr.humans
        mmr.apply_protection(x) ## apply protection if vaccinated
        

        if x.vaccinated 
            @test x.age >= 50 && x.age < 200            
            @test x.protection == 0.95
            v += 1
        else
            @test x.protection == 0.0 
        end
    end
    a1 = length(findall(x -> x.age >= 50 && x.age < 200, mmr.humans))
    #@test isapprox(v/a1, vax_coverage; atol=0.1)

    ## check if a newborn kills the protection level 
    for x in mmr.humans
        if x.vaccinated == true 
            @test x.protection == 0.95
        end
        mmr.newborn(x)
        @test x.protection == 0.0
        @test x.vaccinated == false
        @test x.vaccinetime == 9999
    end

    ## check whether vaccinated, they always stay vaccinated. 
    mmr.init_population()
    mmr.init_vaccination(vax_coverage)
    for t = 1:gridtime
        for x in mmr.humans 
            r = mmr.age_and_death(x)
            vc = mmr.vaccinate(x)
            mmr.apply_protection(x)
            if x.vaccinated 
                @test x.protection = 0.95
            end            
        end
    end
end

function plot_deathdist()
    d = zeros(Int64, 1000)
    for i = 1:1000
        group = rand(agedist)     
        age = rand(agebraks[group])
        d[i] = mmrvaccinedelay.calc_ageofdeath(age)
    end
    return d
end