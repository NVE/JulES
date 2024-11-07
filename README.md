## JulES

JulES is a fundamental energy market simulation model for operational planning, that uses [TuLiPa](https://github.com/NVE/TuLiPa/) as building blocks. JulES clears the power market in small steps with some of the bids generated from deterministic price prognosis models, and stochastic subsystem models for storage valuation. JulES is designed to decompose complex power market problems into subproblems that can be solved faster and with more details in parallel. Each subsystem (e.g. a battery, watercourse or even a hydrogen system) can for example have a storage valuation model and scenario generation that is tailored to the technology they represent.

This is a prototype to test out ideas and inspire the field. Feedback, ideas and contributions are welcome.

### Motivation:
With the transition towards a renewable-based power system, we need models that can represent the new technologies, markets and dynamics. One of [NVE](https://www.nve.no/english/)’s initiatives to improve our understanding of power market modelling, is the ongoing research project called “Power market modelling in Julia”. The goal of this project is to test a new fundamental energy market model for operational planning. The model should be able to simulate the Northern European power market with high temporal resolution, detailed hydropower (or any other technology), uncertainty in weather, and using only open-source software. We want to find out if decomposing the complex power market problem into many smaller subproblems, solving many of them deterministically and with open-source solvers, can give fast and good results. The simulation concept is inspired by how power dispatch is planned in real life, with longer term price prognosis, calculation of individual storage values (e.g. water, battery or gas storage values) with different models for different technologies, and at the end a market clearing algorithm that takes all the details into account.

### JulES

#### Simulation concept

JulES clears the power market in small steps with some of the bids generated from deterministic price prognosis models, and stochastic subsystem models for storage valuation.

The simulation model uses rolling or shrinkable horizons where the underlying models are solved for each time step:
1.	Deterministic price prognosis models generate long-, medium- and short-term prices for different scenarios. The short-term problems have the most details, but they still have a simplified representation of some of the power system elements (i.e. aggregated hydro per price area).
2. The prices for different scenarios are then used in stochastic subsystem models for storage valuation. Each subsystem (e.g., battery or watercourse) can have its own storage valuation model that is tailored to the storage system they represent, for example, short horizons and short-term prices for small batteries and longer horizons and long-term prices for watercourses with multi-year storages. The subsystem models can calculate storage values for each individual storage without the need for advanced end values. This should scale well when we want to solve a complex and detailed power market problem since the smaller detailed subsystems can be solved fast and in parallel. Adding other technologies like hydrogen or gas storage systems should also scale well, since they can be treated as new subsystems that are run in parallel. At the moment, we have two alternative methods for the subsystem storage valuation problem, and it is possible to add more:
   * StochSubsystem: Solve the subsystem problem as a stochastic two-stage LP problem with Benders decomposition.
   * EVPSubsystem: First, solve deterministic LP problems for each scenario for the subsystem, and then use end-values from these problems in a stochastic LP problem with Benders decomposition (same method as StochSubsystem). This combination decomposes the problem even more, allowing for longer horizons, more details, and faster solve times compared to only using StochSubsystem.
3.	Then we clear the power market for one or two days at a time, considering storage values (water and battery values) and end states calculated in models 1 and 2. In the market clearing problem the power system elements can have all their details since the problem horizon is so short.
4. There is also static or dynamic scenario generation at different stages in the JulES algorithm:
   * Reduce the scenarios in the dataset down to a number that should be used throughout the simulation (static / done once).
   * Choose subsystem scenarios from the price static scenarios (dynamic)
        - In the stochastic subsystem models we use scenario modelling to consider uncertainty from all scenarios (30 in the demo) with only a few (7 in the demo). Scenarios can be chosen and weighted with different methods. In the demo we use InflowClusteringMethod which cluster together scenarios with similar total energy inflows in the whole dataset (both level and profile). One scenario from each cluster will represent the others with the weight based on the size of the cluster.
        - At the moment we use the same scenarios for all the subsystems. This works ok since the most important subsystem models are the watercourses in the Nordics. But the scenarios chosen will work worse for the battery systems. In the future we would like to have different scenario modelling for different technologies and systems in different geographical areas. Other alternatives could be to do the scenario modelling based on the price series, or with the residual load (including energy inflow).
5. We have also integrated inflow models into the JulES algorithm, which generate inflow series for each watercourse with state information like snow reservoir levels and weather forecasts. In the future we plan to expand this feature to generate forecasts for demand, wind and solar aswell.

#### Some features we want to highlight
- Water values (and battery storage values) are calculated for each individual storage, without the need for advanced end values.
- The storage valuation model for each subsystem can be customized depending on the technology and geographical location. For example the horizon length and time resolution of the model. We later want to add the possibility to do scenario modelling for each subsystem.
- Since the dataset is stored as time series we can run the model with different time resolutions without having to adapt the dataset.
- Multi-year storage can be considered.
- Scenarios and subsystems are run in parallel.
- Model configuration consists of choosing horizons / temporal resolution and degree of detail for each technology in the dataset, for each type of subproblem.
- Parallel processing, solver warm start, reuse of cuts in multiple time steps, ShrinkableHorizon and scenario generation can be used to make JulES run faster.
- The hydropower is already modelled quite detailed with PQ-curves, environmental constraints, ramping restrictions and head dependency.
- Scenarios are phased in from the main/simulation scenario since uncertainty is lower close to the decision time.
- The price prognosis, subsystem and clearing problem uses the same horizon, a horizon with a high time resolution in the start and a decreasing granularity. This simplifies the transfer of data between the different models. As an example the clearing problem will only use the first days of the horizon, the stochastic subsystem problem the first two months, the end value problem the first year, and the price prognosis problem will use the full horizon of maybe 5 years. 

#### TuLiPa
TODO: Benefits of [TuLiPa](https://github.com/NVE/TuLiPa/)

#### Get an overview of JulES:
- src/run_serial.jl - Code for running JulES simulation
- src/prob_ppp.jl – Code for price prognosis problems
- src/prob_evp.jl – Code for end value problems
- src/prob_stoch.jl – Code for stochastic subsystem problems
- src/prob_cp.jl - Code for market clearing problem
- src/scenariomodelling.jl - Code for scenario modelling
- src/ifm.jl - Code for inflow model
- src/local_db.jl - Code for in-memory process-local database
- src/util.jl - Various useful functions

#### See also demos (&#x2714; = open data so you can run it yourself):
- NB!!! For this release we have not prioritized the demos, so they have no comments.
- demos/Demo JulES as a solar with battery model &#x2714;
- demos/Demo JulES as a single watercourse model &#x2714;
- demos/Demo JulES as a long-term series simulation model
- demos/Demo JulES as a medium-term parallel prognosis model

#### Results:
We are at the moment testing JulES as a long-term series simulation model, a medium term parallel prognosis model, a single watercourse model, and as a solar with battery model. We will publish results from this work at a later stage.

### Setup
*  Install julia version 1.9.2:
https://julialang.org/downloads/
* Clone repository
```console
git clone https://github.com/NVE/JulES.git
``` 

Enter the folder with project.toml using a terminal window and run Julia and these commands:

With julia prompt change to Julias Pkg mode using ] and enter

```console
Julia> ]
```
With pkg prompt and while being inside the project folder, activate the project
```console
(@v1.9) pkg> activate .
```

Prompts shows the project is activated, then installs the libraries needed with instantiate
```console
(JulES) pkg> instantiate
```

To start running the demos, run jupyter notebook from the terminal while being inside the JulES project folder (make sure the kernel is julia 1.9.2)
```console
jupyter notebook 
```

If jupyter is not installed then it can be installed using IJulia
```console
Julia> using IJulia
Julia> IJulia.notebook()
```

### Contact:
Julien Cabrol: jgrc@nve.no

Harald Endresen Haukeli: haen@nve.no

### Licensing:
Copyright 2023 The Norwegian Water Resources and Energy Directorate, and contributors.

JulES is free software; you can redistribute it and/or modify it under the terms of the MIT license. See [COPYING](https://github.com/NVE/JulES/blob/master/COPYING) for the license text.
