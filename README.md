## JulES

JulES is a fundamental energy market simulation model for operational planning, that uses [TuLiPa](https://github.com/NVE/TuLiPa/) as building blocks. JulES clears the power market in small steps with some of the bids generated from deterministic price prognosis models, and stochastic subsystem models for storage valuation. Jules is designed to decompose complex power market problems into subproblems that can be solved faster and with more details in parallel. Each subsystem (e.g. a battery, watercourse or even a hydrogen system) can for example have a storage valuation model and scenario generation that is tailored to the technology they represent.

This is a prototype to test out ideas and inspire the field. Feedback, ideas and contributions are welcome.

### Motivation:
With the transition towards a renewable-based power system, we need models that can represent the new technologies, markets and dynamics. One of [NVE](https://www.nve.no/english/)’s initiatives to improve our understanding of power market modelling, is the ongoing research project called “Power market modelling in Julia”. The goal of this project is to test a new fundamental energy market model for operational planning. The model should be able to simulate the Northern European power market with high temporal resolution, detailed hydropower (or any other technology), uncertainty in weather, and using only open-source software. We want to find out if decomposing the complex power market problem into many smaller subproblems, solving many of them deterministically and with open-source solvers, can give fast and good results. The simulation concept is inspired by how power dispatch is planned in real life, with longer term price prognosis, calculation of individual storage values (e.g. water, battery or gas storage values) with different models for different technologies, and at the end a market clearing algorithm that takes all the details into account.

### JulES

#### Simulation concept

JulES clears the power market in small steps with some of the bids generated from deterministic price prognosis models, and stochastic subsystem models for storage valuation.

The simulation model uses rolling or shrinkable horizons where the underlying models are solved for each time step:
1.	Deterministic price prognosis models generate long-, medium- and short-term prices for different scenarios. The short-term problems have the most details, but they still have a simplified representation of some of the power system elements (i.e. aggregated hydro per price area).
2.	The prices for different scenarios are then used in stochastic subsystem models for storage valuation. Each subsystem (e.g. battery or watercourse) can have its own storage valuation model that is tailored to the storage system they represent, for example short horizons and short-term prices for small batteries and longer horizons and long-term prices for watercourses with multi-year storages. The subsystem models can calculate storage values for each individual storage without the need for advanced end values. This should scale well when we want to solve a complex and detailed power market problem since the smaller detailed subsystems can be solved fast and in parallel. Adding other technologies like hydrogen or gas storage systems should also scale well, since they can be treated as new subsystems that are run in parallel.
3.	Then we clear the power market for one or two days at a time, considering storage values (water and battery values) and end states calculated in models 1 and 2. In the market clearing problem the power system elements can have all their details since the problem horizon is so short.
4. There is also static or dynamic scenario generation at different stages in the JulES algorithm:
   * Reduce the scenarios in the dataset down to a number that should be used in the simulation (static / done once).
   * Choose price prognosis scenarios from the simulation scenarios (dynamic)
   * Choose subsystem scenarios from the price prognosis scenarios (dynamic)
        - In the stochastic subsystem models we use scenario modelling to consider uncertainty from all scenarios (30 in the demo) with only a few (7 in the demo). Scenarios can be chosen and weighted with different methods. In the demo we use InflowClusteringMethod which cluster together scenarios with similar total energy inflows in the whole dataset (both level and profile). One scenario from each cluster will represent the others with the weight based on the size of the cluster.
        - At the moment we use the same scenarios for all the subsystems. This works ok since the most important subsystem models are the watercourses in the Nordics. But the scenarios chosen will work worse for the battery systems. In the future we would like to have different scenario modelling for different technologies and systems in different geographical areas. Other alternatives could be to do the scenario modelling based on the price series, or with the residual load (including energy inflow.)
        - We also want to implement more complex scenario generation. This modular design of JulES gives us the flexibility to in the future generate inflow series for each watercourse with state information like snow reservoir levels and weather forecasts.

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

#### TuLiPa
TODO: Benefits of [TuLiPa](https://github.com/NVE/TuLiPa/)

#### Get an overview of JulES:
- src/prognosis.jl – Code for price prognosis models
- src/stochastic.jl – Code for stochastic sub system models
- src/clearing.jl - Code for market clearing problem
- src/scenariomodelling.jl - Code for scenario modelling
- src/util.jl - Various useful functions

#### See also demos (&#x2714; = open data so you can run it yourself):
- NB!!! For this release we have not prioritized the demos, so they have no comments.
- demos/Demo la.ipynb
- demos/Demo prognosis.ipynb
- demos/Demo solar battery.ipynb :heavy_check_mark:
- demos/Demo watercourse.ipynb :heavy_check_mark:

#### Results:
We are at the moment testing JulES as a medium term parallel prognosis model. The results can then be compared to our other prognosis model (EMPS) and quickly be compared to the real market and historical data. We will publish results from this work at a later stage.

We have also tested JulES on a simplified (most importantly simple transmission modelling) version of the Nortwestern European dataset NVE uses for its long-term analyses. This dataset contains:
  - 32 price areas (9 exogen)
  - 73 transmission lines (19 with ramping)
  - 162 demands
  - 88 batteries (normal, and representing demand response and V2G)
  - 294 thermal plants (228 with start up costs)
  - 100 wind and solar plants (aggregated)
  - 1482 hydropower modules
    - 965 with production
    - 43 with pumps (includes PHS)
    - 998 with reservoirs
    - 788 restrictions (environmental, reservoir curves and ramping)
    - 90 PQ-curves (mostly Sweden)
    - Metadata for head dependency (nominal head, outlet level and reservoir head to filling curves) for some plants and pumps
      
In the first phase of JulES we have prioritized implementing functionality and test that they give the intended results, and compared results against other long-term models (TheMA and Samnett, which give a sufficient basis for comparison) on a big dataset. The price formation, reservoir operation and runtime of JulES seems promising, but there remains a lot of testing to see all the impacts of the concept implementation, model parameters, scenario modelling and head dependency on water values, prices, reservoir operation and runtime. The testing has been restricted by the runtime of the model, which we will have to improve to test different parts of the model more efficiently.
  
Throughout the testing we have achieved the wanted price volatility in the thermal dominated part of the dataset (Western Europe). On the other hand, the Nordics have had very flat prices due to too much flexibility in the hydropower modelling. We have therefore added hydropower production ramping restrictions in the market clearing in an attempt to reduce the flexibility of the run-of-river hydropower plants. This results in much more price volatiliy, but at a big computational cost. Ramping restrictions on transmission lines has the same effect.

The runtime of the model is promising considering the big dataset, and the possible optimization we have in mind:
* We have already implemented and tested the commercial solver CPLEX (we now use HiGHS), which especially improves runtime for the market clearing problem.
* We have implemented ShrinkableHorizon, which reduces update and solve times of the price prognosis problems.
* We are actively working on testing JulES on more and faster processor cores in parallel. In the branch **redesign_distributed** we are changing the structure of JulES to allow more flexibility in the parallelization, and reduce the memory usage of the model. 
* We are not done exploring the model settings, for example duration of horizons and time resolutions in the different problem. For example we could clear the market for 24 hours at a time instead of 48 hours like now, which could reduce the computational time depending on if the market clearing has a higher runtime than the stochastic subsystem and price prognosis models (this is the case when we have a very detailed market clearing).
* We could try different configurations of ramping restrictions, and test if time delays in watercourses can achieve the same effects at a lower computation cost. Considering unavailability of hydropower or reserve market obligations, should also decrease the flexibility of the hydropower system. Detailed transmission system modelling should also be implemented in the future.
* We could run the model only for the Nordics, which would reduce the size of the dataset substantially and give results that are comparable to other models we use.

The price levels in the Nordics are higher than in other models. This is partly due to the high flexibility in the hydropower modelling, which gives stable high prices and not many zero-prices. Another reason is that the stochastic subsystem models could need some improvements, for example longer and more detailed horizons. This should give water values that gives better long term signals. The water value calculations are also being improved in the branch **redesign_distributed**. More load shedding can also be a contributor to higher prices, but this can be prevented with scenario modelling, head dependencies and flexible demand. 

We have also seen the effects of scenario modelling and head dependencies. Scenario modelling can be used to reduce the runtime and adjust the risk taking, which gives more realistic reservoir operation and avoid the extremes of flooding and load shedding. We are also working on expanding the scenario modelling with inflow models that account for the snow level in the current state. This should give a more realistic outcome space for the inflow scenarios. Head dependency can be used to get a more realistic reservoir operation and higher production, and also gives lower risk of load shedding. These have to be tested further. 

We are also very happy with the modelling choice of modularity, using time-series datasets and using Julia. This has made TuLiPa and JulES very pleasant to work with, as they provide a great deal of flexibility in adding complex functionality without having to make extensive changes to the existing code. Additionally, the models can be run with different methods and time resolutions without adaptations of the dataset. These design choices contribute to the model's suitability for further development and modeling the future power system when new modeling requirements arise.

However, the project's codebase needs to be professionalized with better structure to make the model more user-friendly, allowing not only developers to run the model. Unit testing is also important to ensure that the model functions as intended. So far, we have been working on the model concept alone, so it will be crucial to involve analysts who will use the models and developers outside of NVE who can contribute to further developing the concept.

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
Julia> Julia.notebook()
```

### Contact:
Julien Cabrol: jgrc@nve.no

Harald Endresen Haukeli: haen@nve.no

### Licensing:
Copyright 2023 The Norwegian Water Resources and Energy Directorate, and contributors.

JulES is free software; you can redistribute it and/or modify it under the terms of the MIT license. See [COPYING](https://github.com/NVE/JulES/blob/master/COPYING) for the license text.
