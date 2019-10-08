module mmrvaccinedelay

using Parameters      ## with julia 1.1 this is now built in.
using ProgressMeter   ## can now handle parallel with progress_pmap
#using PmapProgressMeter
#using DataFrames
using Distributions
using StatsBase
using StaticArrays
using Random

@enum HEALTH UNDEF = 0 SUSC = 1 INF = 2 REC = 3 VAC = 4 

mutable struct Human ## mutable structs are stored on the heap 
    idx::Int64
    health::HEALTH
    swap::HEALTH #do we need  this? We can do a sequential atualization
    
    age::Int64     
    group::Int64      

    Human() = new()
end

#global parameters 
const agedist =  Categorical(@SVector [0.053, 0.055, 0.052, 0.056, 0.067, 0.07, 0.07, 0.068, 0.064, 0.066, 0.072, 0.073, 0.064, 0.054, 0.042, 0.029, 0.045])
const agebraks = @SVector [0:4, 5:9, 10:14, 15:19, 20:24, 25:29, 30:34, 35:39, 40:44, 45:49, 50:54, 55:59, 60:64, 65:69, 70:74, 75:79, 80:85]

const humans = Array{Human}(undef, 10000)
export HEALTH, humans


function init_human(x::Human, idx)
   init_human(x)
   x.idx = idx         
end

function init_human(x::Human)
    x.idx = 0
    x.health = SUSC
    x.swap = UNDEF 
    g = rand(agedist) ## get an agegroup
    x.age = rand(agebraks[g])
    x.group = min(g, 15)         
end
export init_human

function replace_human(h::Human)
    oldid =  h.idx
    oldgrp = h.contact_group
    init_human(h)
    h.age = 0
    h.idx = oldid
end
export replace_human

function init_population()    
    @inbounds for i = 1:10000
        humans[i] = Human()   ## create an empty human
        init_human(humans[i], i) ## initialize the human
    end
end
export init_population


greet() = print("Hello World!")

end # module
