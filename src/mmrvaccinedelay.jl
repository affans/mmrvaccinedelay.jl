## to do, 
## convert contact matrix to PDF/Categoricals.

module mmrvaccinedelay

# 4 major components to add
# weekly unit step, for vaccination delay purposes. 
# contact patterns. run full test suite
# seasonal beta, use some sort of function with week input
# age groups, for vaccination and calibration purposes

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
    swap::HEALTH     
    age::Int64          # in weeks 
    ageofdeath::Int64   # in weeks
    group::Int64     
    vaccinated::Bool
    vaccinetime::Int64  # in weeks 
    Human() = new()
end

#global parameters 
# agebraks not really used
const agedist =  Categorical(@SVector [0.053, 0.055, 0.052, 0.056, 0.067, 0.07, 0.07, 0.068, 0.064, 0.066, 0.072, 0.073, 0.064, 0.054, 0.116])
const agebraks_old = @SVector [0:4, 5:9, 10:14, 15:19, 20:24, 25:29, 30:34, 35:39, 40:44, 45:49, 50:54, 55:59, 60:64, 65:69, 70:85]
const agebraks = @SVector [0:200, 201:450, 451:700, 701:950, 951:1200, 1201:1450, 1451:1700, 1701:1950, 1951:2200, 2201:2450, 2451:2700, 2701:2950, 2951:3200, 3201:3450, 3451:4250]

const humans = Array{Human}(undef, 10000)
export HEALTH, humans

const NB = negative_binomials()
const CM = contact_matrix()

Base.show(io::IO, ::MIME"text/plain", z::Human) = dump(z)

function main(simnumber)
    # main entry of the simulation
    Random.seed!(simnumber)
    init_population()

    insert_infected(2000)



    maxtime = 50*50
    
    ## data collection variables = number of elements is the time units. 
    ## so the vector collects number of latent/symp/asymp at time t. 
    ## it does not collect the initial latent case.
    vacc_ctr = zeros(Int64, maxtime )   
    symp_ctr = zeros(Int64, maxtime)   

    ## contact matrices

    ## these matrices are used to calculate contact patterns
    fcm = zeros(Int64, 15, 15)   ## how many times did susc/sick contact group i meet contact group j meet but failed to infect.
    cmg = zeros(Int64, 15, 15)   ## how many times did contact group i meet with contact group j
    
    ## we setup an array for who belongs to each group
    _contacts = zeros(Int64, 10000, 50*50)
    
    infected = findall(x -> x.health == INF, humans)
    oldc = _contacts
    for i = 1:maxtime
        agm = [humans[i].group for i = 1:10000]

        for j in 1:10000
            g, c = age_and_vaccinate(humans[j])  ## increase age of each agent. vaccinate if neccessary            
            if g 
                #println("g changed")
                edit_contacts(j, _contacts)
            end
        end
        # for i in infected
        #     contact_dynamic2(humans[i], NB, CM, agm)
        # end
    end    
    return oldc, _contacts
end
export main

function setup_contacts(cts)
    for i = 1:10000
        ig = humans[i].group                
        cts[i, :] = rand(NB[ig], 50*50) ## get the number of contacts
    end
    return cts
end
export setup_contacts

@inline function edit_contacts(rowid, cts)
    cts[rowid, :] = rand(NB[humans[rowid].group ], 50*50)
end
export edit_contacts

function insert_infected(cnt) 
    h = rand(1:10000, cnt)
    for i in h 
        humans[i].health = INF
    end
end
export insert_infected

function init_human(x::Human, idx = 0)
    x.idx = idx
    x.health = SUSC
    x.swap = UNDEF 
    x.group = rand(agedist)     
    x.age = rand(agebraks[x.group])
    # age of death in group 15 
    x.ageofdeath = rand(3500:3750)   ## could happen that x.ageofdeath < x.age in which case they will die right away
    x.vaccinated = false
    x.vaccinetime = 9999 
end
export init_human

function replace_human(h::Human)
    oldid =  h.idx
    init_human(h)
    h.age = 0
    h.group = 1
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

function get_group(age)
    # input age: in weeks, output group from 1:15 
    r = findfirst(brak -> age in brak, agebraks) 
    if r === nothing 
        println("$age")
        error("get_group threw an error")
    end
    return r 
end
export get_group

function age_and_vaccinate(x::Human)
    grp_changed = false
    cnt_vaccinated = 0
    grp = x.group
    newage = x.age + 1

    if newage > x.ageofdeath
        replace_human(x)
        ls = true
    else 
        x.age = newage
        x.group = get_group(newage) 
    end
    # group has changed
    if x.group != grp 
        grp_changed = true
    end

    if newage == x.vaccinetime
        x.vaccinated == true
        cnt_vaccinated += 1
    end
   
    return grp_changed, cnt_vaccinated
end
export age_and_vaccinate


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
export contact_matrix

function negative_binomials()
    ##Age group's mean
    AgeMean = Vector{Float64}(undef, 15)
    AgeSD = Vector{Float64}(undef, 15)

    AgeMean = [10.21, 14.81, 18.22, 17.58, 13.57, 13.57, 14.14, 14.14, 13.83, 13.83, 12.3, 12.3, 9.21, 9.21, 6.89]
    AgeSD = [7.65, 10.09, 12.27, 12.03, 10.6, 10.6, 10.15, 10.15, 10.86, 10.86, 10.23, 10.23, 7.96, 7.96, 5.83]

    nbinoms = Vector{NegativeBinomial{Float64}}(undef, 15)
    for i = 1:15
        p = 1 - (AgeSD[i]^2-AgeMean[i])/(AgeSD[i]^2)
        r = AgeMean[i]^2/(AgeSD[i]^2-AgeMean[i])
        nbinoms[i] =  NegativeBinomial(r, p)
    end

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
export negative_binomials

function contact_dynamic2(x, NB, CM, agm)
    ## right away check if x is infected, otherwise no use.
            
    ig = x.group            
    _contacts = rand(NB[ig]) ## get the number of contacts
    ds = Categorical(CM[ig])
    _groups = rand(ds, _contacts) ## from which group.

    u=unique(_groups)
    d=Dict([(i,count(x->x==i,_groups)) for i in u])
    for (k, v) in d
        possagents = findall(x -> x == k, agm)
        meetagents = rand(possagents, v) 
        for j in meetagents
            if humans[j].health == SUSC 
                if rand() < 0.5 
                    #println("beta not implemented yet")
                    humans[j].health = INF
                end
            end 
        end
    end    
end
export contact_dynamic2



greet() = print("Hello World edited!")
export greet

end # module
