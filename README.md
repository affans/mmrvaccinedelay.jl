# MMR Vaccine Delay Model

## Evaluating the effects of delaying the mmr vaccine.
----
## Reproducibility of the model

### How to download
The model is self-contained (non-registered) Julia package. To use the model, add the package by typing `] add https://github.com/affans/mmrvaccinedelay.jl` in the Julia REPL.

Once the package is downloaded, a single *realization* of the model can be run by invoking the `main` function. This function accepts a few arguments (settings of the model). Description of the function arguments can be found by typing `? mmrvaccinedelay.main` in the Julia REPL.  

Since the model is stochastic, many realizations are required to get an accurate picture of the results. This essentially means running the `main` function repeatedly, saving the results for each replicate. This can be done very easily using Julia's Distributed library. Simply `addprocs` or using `ClusterManagers` to connect to a cluster to launch `n` number of parallel workers. Then use `pmap` to run the main function on each worker, saving the results to be processed later. 

Example of this multitasking can be seen in `initcluster.jl` and `process.jl` files. 

For the advanced users, you may `dev` the package instead (see "Program Structure" for an overview of the model). In addition, activating the correct package environment allows for `] test` for unit-testing the model. 

### Program Structure
This section is mainly for me but also anyone else who is interested in modifying the code. It explains the main functions of the model. 

The main program begins in `main()`. 
After the parameters are set up, three main functions are run: 

1. `init_population` initializes the population by calling `apply_agedistribution` which distributes the age (and calculates age of death), and `apply_protection`. 

2. `init_herdimmunity` inserts a portion of `REC` individuals who have 100% protection. 
3.  `insert_infected` and inserts some number (possibly zero) random infected individual at the start of simulations.

After the model is initialized, the main time loop starts. The time variable `t` is checked against a pre-specified time for starting vaccine, at which `P.vaccine_onoff = true` for new borns and applies vaccination of existing 1 - 4 years olds (with some pre-specified coverage level).

Next, if the time variable `t` is equal to a pre-specified time for an outbreak, the `insert_infected` function is used to introduce 1 or more infected individuals which is expected to cause an outbreak. 

Next, the model loops over every infected individual and checks for transmission using  `contact_dynamics()`. If disease is transfer, a `swap` is set, and later the health status of the individual is updated using the `update_swaps` function. 

Next, the model runs the `age_and_death()` function which increases the age of an individual. If the age reaches the age of death, the person dies and a newborn is put in their place. 

Next, the model runs the `vaccinate()` function for each individual. If it is the person's time to vaccinate, the individual will have `x.vaccinated` property to true. 

Next, `apply_protection()` is run to calculate the correct protection level of all individuals. 

The loop ends and all the counters are returned from the function for analysis. 








