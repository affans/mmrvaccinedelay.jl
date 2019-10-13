using mmrvaccinedelay
using Test

const mmr = mmrvaccinedelay
const gridsize = length(mmr.humans)
const gridtime = 50*50
const vax_coverage = 0.8
@testset "INIT" begin
    mmr.init_population()
    x = mmr.humans[1]
    @test x.idx == 1
    @test x.vaccinetime == 9999
    @test x.vaccinated == false
    @test length(findall(x -> x == undef, mmr.humans)) == 0 ## check for proper initialization
    for x in mmr.humans
        @test x.ageofdeath == 4250
        @test x.swap == mmr.UNDEF
        @test x.health == mmr.SUSC
    end
    
    ## test if the init_infected function really works.
    mmr.init_population()
    mmr.init_infected(1)
    a = findall(x -> x.health == mmr.INF, mmr.humans)
    @test length(a) == 1

    ## test if the init_infected dosn't sample the sample person twice.
    mmr.init_population()
    mmr.init_infected(100)
    a = findall(x -> x.health == mmr.INF, mmr.humans)
    @test length(a) == 100


    mmr.init_population()
    for x in mmr.humans 
        @test x.protection == 0.0
        @test x.health == mmr.SUSC
    end
    mmr.init_herdimmunity()
    a1 = findall(x -> x.health == mmr.SUSC && x.age < 26, mmr.humans)
    a2 = findall(x -> x.health == mmr.SUSC && x.age >= 26 && x.age < 50, mmr.humans)
    a3 = findall(x -> x.health == mmr.REC, mmr.humans)
    
    for i in a1 
        @test mmr.humans[i].protection == 1.0
        @test mmr.humans[i].vaccinated == false ## extra test
    end
    for i in a2
        @test mmr.humans[i].protection == 0.5
        @test mmr.humans[i].vaccinated == false ## extra test
    end
    for i in a3
        @test mmr.humans[i].protection == 1.0
        @test mmr.humans[i].vaccinated == false ## extra test
    end
    p = findall(x -> x.age >= 200, mmr.humans)
    @test isapprox(length(a3)/length(p), 0.80; atol=0.1)
    
end

@testset "AGE" begin
    mmr.init_population()    
    ## check if age_and_vaccinate function (age only)
    ## checks if age is incremented properly
    ## and verifies dead people are replaced properly. 
    a1 = [mmr.humans[i].age for i = 1:gridsize]
    ls = zeros(Bool, gridsize) #contains whether they have left the system
    for i = 1:gridsize
        r = mmr.age_and_vaccinate(mmr.humans[i])
        ls[i] = r[2] # r[2] contains whether they have left the system
    end
    a2 = [mmr.humans[i].age for i = 1:gridsize]
    for i = 1:gridsize
        @test (a2[i] == a1[i] + 1) || a2[i] == 0
    end
    l1 = length(findall(x -> x == 0, a2))
    l2 = length(findall(x -> x == true, ls))
    @test l1 == l2

    ## check whether group changes (may as well loop through all humans, all groups)
    mmr.init_population()
    for x in mmr.humans
    #x = mmr.humans[1]
    for grp = 1:15 
        # set their age to the end of the group and see if age and group changes accordingly
        x.age = mmr.agebraks[grp][end]
        x.group = mmr.get_group(x.age) 
        
        @test x.group == grp
        @test x.ageofdeath == 4250 ## see if this dosn't change
        oldage = x.age
        oldgrp = x.group
        r = mmr.age_and_vaccinate(x)
        @test r[1] == true  ## test if the age_and_vaccinate returns true for group change
     
        if grp < 15  ## if group is less than 15 .. they won't be reset
            @test x.age == oldage + 1
            @test x.group == oldgrp + 1
            @test x.group == mmr.get_group(x.age) ## this tests the function
        else # since group was 15 and their age was set to the last age
            @test x.age == 0
            @test x.group == 1
            @test x.group == mmr.get_group(x.age)
            @test r[2] == true ## r[2] holds the true/value whether person has left.
        end

        # set their age to the START of the group and see if age and group changes accordingly 
        # group should not change, no one should die here either
        x.age = mmr.agebraks[grp][1]
        x.group = mmr.get_group(x.age) 
        oldage = x.age
        r = mmr.age_and_vaccinate(x)
        @test r[1] == false
        @test r[2] == false
        @test x.group == grp 
        @test (x.age == oldage + 1)
    end
    end

    ## this is not a test. basically runs the model for gridtime and 
    ## stores who left the system. this needs to be plotted or something
    mmr.init_population()
    wd = zeros(Int64, gridtime)
    for t = 1:gridtime
        for x in mmr.humans 
            r = mmr.age_and_vaccinate(x)
            if r[2] == true 
                wd[t] += 1
            end
        end
    end
    #UnicodePlots.scatterplot(wd)

    ## we test that if everywhere was born then within 50 years they should not leave the system
    ## secondary test that ID of the person should not be changing. the .idx proeprty of human is important
    mmr.init_population()
    wd = zeros(Int64, gridtime)
    for x in mmr.humans 
        x.age = 0 
        x.group = 1
    end
    for t = 1:gridtime
        for x in mmr.humans 
            oldid = x.idx
            r = mmr.age_and_vaccinate(x)
            if r[2] == true
                @test x.age == 0 
                @test x.group == 1
                @test oldid == x.idx 
                wd[t] += 1
            end
        end
    end
    @test sum(wd) == 0
   
    # test agm vector. this is used in the contact patterns. 
    mmr.init_population()
    x = mmr.humans[1]
    x.age = 0
    x.group = 1
    _agm = Vector{Vector{Int64}}(undef, 15)
    mmr.setup_agm(_agm)
    @test 1 ∈ _agm[1] ## test the manual human above. 
    for x in mmr.humans  # test everyone
        @test x.idx ∈ _agm[x.group]         
    end
    x.age = 201
    x.group = mmr.get_group(x.age)
    mmr.setup_agm(_agm)
    @test 1 ∉ _agm[1] # human changed group so should not in the first bin
    @test 1 ∈ _agm[2]
end

@testset "SWAPS" begin

    ## This tests whether the person's swaps are set properly.
    ## secondary tests makes sure if person isn't infected, they can't transmit the disease even with beta = 1
    mmr.init_population()
    x = mmr.humans[1]            # find a human
    x.swap = mmr.INF             # set SWAP to infected
    mmr.update_swaps()           # update swap
    @test x.swap == mmr.UNDEF    # tests whether swap happened
    @test x.health == mmr.INF
    # prep for the loop
    _agm = Vector{Vector{Int64}}(undef, 15)
    mmr.setup_agm(_agm)
    # go through one week and see what happens after 
    for i in 1:gridsize
        cntinf = mmr.contact_dynamic2(mmr.humans[i], 1.0, _agm)
        # we only made the first human infected, so the number of secondary infections should be zero for everyone else
        # @test cntinf > 0 ## this could fail.... if number of contacts is zero      
    end 
    @test x.health == mmr.INF    # human[1] is still infected 
    mmr.update_swaps()           # this should move their health to REC
    @test x.health == mmr.REC    # test whether it does
    @test x.swap == mmr.UNDEF
    # test the population to see only one recovered
    @test 1 == length(findall(x -> x.health == mmr.REC, mmr.humans))
    # test the population to find all other infecteds. 
    @test 0 == length(findall(x -> x.swap != mmr.UNDEF, mmr.humans)) # no one should have a swap set
    @test length(findall(x -> x.health == mmr.INF, mmr.humans)) >= 0 # no one should have a swap set
end

@testset "VAX" begin
    ## test that all vaccine time is above 52 weeks
    mmr.init_population()
    mmr.init_vaccination()
    v = 0
    for x in mmr.humans
        if x.vaccinated 
            @test x.age >= 50 && x.age < 200
            @test x.health == mmr.SUSC
            v += 1
        end
    end
    a1 = length(findall(x -> x.age >= 50 && x.age < 200, mmr.humans))
    @test isapprox(v/a1, vax_coverage; atol=0.1)

    mmr.init_population()
    mmr.init_herdimmunity()
    mmr.init_vaccination()
    for t in 1:gridtime
        for x in mmr.humans
            mmr.age_and_vaccinate(x)
        end
    end
    
end
