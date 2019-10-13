## this file just sets up parallel processing AND imports packages
## you can modify this file to run `addprocs` based  
## on available computational structure and capacity.

using Distributed
using Random
using DelimitedFiles
using Distributions
using Base.Filesystem
using DataFrames
using CSV
using ClusterManagers
using Query
using Statistics
using UnicodePlots
using Dates

using mmrvaccinedelay

## multi core. 
##addprocs(4; exeflags="--project=.")
## a = SlurmManager(200, N=17)
#addprocs(SlurmManager(10*32), N=32) ## this throws an error because N = 32. Create PR. 
addprocs(SlurmManager(512), N=16) 

@everywhere using mmrvaccinedelay
@everywhere using ProgressMeter