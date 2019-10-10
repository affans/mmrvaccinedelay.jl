using mmrvaccinedelay
using Test

const mmr = mmrvaccinedelay
const l = length(mmr.humans)
@testset "initialization" begin
    mmr.init_population()
    x = humans[1]
    @test x.idx == 1
    @test x.vaccinetime == 9999
    @test x.vaccinated == false
    @test length(findall(x -> x == undef, humans)) == 0 ## check for proper initialization
end

@testset "age&vaccinate" begin
    
    mmr.init_population()    
    ## check if age_and_vaccinate function (age only)
    ## checks if age is incremented properly
    ## and verifies dead people are replaced properly. 
    a1 = [humans[i].age for i = 1:10000]
    ls = zeros(Bool, l)
    for i = 1:l
        ls[i] = mmr.age_and_vaccinate(humans[i])
    end
    a2 = [humans[i].age for i = 1:10000]
    for i = 1:l
        @test (a2[i] == a1[i] + 1) || a2[i] == 0
    end
    l1 = length(findall(x -> x == 0, a2))
    l2 = length(findall(x -> x == true, ls))
    @test l1 == ls

    ## to do
    ## check whether groups are changing as age++
    ## repeat the loop for 1000 times and see how age changes (store it into a matrix)
    ## check whether agm contains the right groups. 
    ## check whether the sum of all IDs is the same after 1000 iterations. the .idx proeprty of human is important
end
