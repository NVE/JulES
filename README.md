## JulES

JulES is an energy market simulation model for long term planning that uses [TuLiPa](https://github.com/NVE/TuLiPa/) as building blocks. JulES clears the power market in small steps with some of the bids generated from deterministic price prognosis models, and stochastic subsystem models for storage valuation. Jules is designed to decompose complex power market problems into subproblems that can be solved faster and with more details in parallel. Each subsystem (e.g. a battery, watercourse or even a hydrogen system) can for example have a storage valuation model that is tailored to the technology they represent.

This is a prototype to test out ideas and inspire the field. Feedback, ideas and contributions are welcome.

### Motivation:
With the transition towards a renewable-based power system, we need models that can represent the new technologies, markets and dynamics. One of [NVE](https://www.nve.no/english/)’s initiatives to improve our understanding of power market modelling, is the ongoing research project called “Power market modelling in Julia”. The goal of this project is to make an algorithm for simulating the Northern European power market with high temporal resolution, detailed hydropower (or any other technology), uncertainty in weather, and using only open-source software. We want to find out if decomposing the complex power market problem into many smaller subproblems, solving many of them deterministically and with open-source solvers, can give fast and good results. The simulation concept is inspired by how power dispatch is planned in real life, with longer term price prognosis, calculation of individual storage values (e.g. water, battery or gas storage values) with different models for different technologies, and at the end a market clearing algorithm that takes all the details into account. 

### JulES

#### Simulation concept

JulES clears the power market in small steps with some of the bids generated from deterministic price prognosis models, and stochastic subsystem models for storage valuation.

The simulation model uses a rolling horizon approach where the underlying models are solved for each time step:
1.	Deterministic price prognosis models generate long-, medium- and short-term prices for different scenarios. The short-term problems have the most details, but they still have a simplified representation of some of the power system elements (i.e. aggregated hydro per price area).
2.	The prices for different scenarios are then used in stochastic subsystem models for storage valuation. Each subsystem (e.g. battery or watercourse) can have its own storage valuation model that is tailored to the storage system they represent, for example short horizons and short-term prices for small batteries and longer horizons and long-term prices for watercourses with multi-year storages. The subsystem models can calculate storage values for each individual storage without the need for advanced end values. This should scale well when we want to solve a complex and detailed power market problem since the smaller detailed subsystems can be solved fast and in parallel. Adding other technologies like hydrogen or gas storage systems should also scale well, since they can be treated as new subsystems that are run in parallel. 
3.	Then we clear the power market for one or two days at a time, considering storage values (water and battery values) and end states calculated in models 1 and 2. In the market clearing problem the power system elements can have all their details since the problem horizon is so short.

#### Some features we want to highlight
- Scenarios and subsystems are run in parallel.
- The storage valuation model for each subsystem can be customized depending on the technology and geographical location. For example the horizon length, time resolution or scenario generation of the model.
- Water values (and battery storage values) are calculated for each individual storage, without the need for advanced end values.
- Multi-year storage can be considered.
- Model configuration consists of choosing horizons / temporal resolution and degree of detail for each technology in the dataset, for each type of subproblem.
- Parallel processing, solver warm start, reuse of cuts in multiple time steps, and scenario generation can be used to make JulES run faster.
- We have implemented scenario reduction and an interface to allow different methods. We later want to implement scenario generation.
- The hydropower is already modelled quite detailed with PQ-curves, environmental constraints and head dependency.

#### TuLiPa
TODO: Benefits of [TuLiPa](https://github.com/NVE/TuLiPa/)

#### Get an overview of JulES:
- src/prognosis.jl – Code for price prognosis models
- src/stochastic.jl – Code for stochastic sub system models
- src/clearing.jl - Code for market clearing problem
- src/scenarioreduction.jl - Code for scenario reduction
- src/util.jl - Various useful functions

#### See also demos:
- demos/Prototype JulES v1

#### Results:
- We have tested JulES on a simplified (most importantly simple transmission modelling) version of the Nortwestern European dataset NVE uses for its analyses . This dataset contains:
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
    - Metadata for head dependency (nominal head, outlet level and reservoir curves) for some plants and pumps
- Throughout the testing we have achieved the wanted price volatility in the thermal dominated part of the dataset (Western Europe). On the other hand, the Nordics have had very flat prices due to too much flexibility in the hydropower system. In this demo we have therefore added hydropower production ramping restrictions in an attempt to reduce the flexibility of the run-of-river hydropower plants. This results in much more price volatiliy, but at a big computational cost. 
- (not updated!!!) Ramping restrictions on transmission lines has the same effect. With ramping restrictions (hydro and transmission) the computational time in the demo is around 10 seconds per simulated day. Without ramping in the market clearing it is down to around 4.5 seconds per day. This is promising considering the big dataset, and the list of possible optimization we have in mind. It is also always possible to clear the market for 24 hours instead of 48 hours like now, which would reduce the computational time substantially. It is also interesting what these computational times would be with a commercial solver (now we use HiGHS).
- We will try different configurations of ramping restrictions, and also test if time delays can achieve the same effects at a lower computation cost. Considering unavailability of hydropower or reserve market obligations, should also decrease the flexibility of the hydropower system.

#### Possible improvements to JulES:
See file "Possible improvements to JulES"

### Contact:
Julien Cabrol: jgrc@nve.no

Harald Endresen:

### Licensing:
Copyright 2023 The Norwegian Water Resources and Energy Directorate, and contributors.

Specific license not decided, but the code will be open-source.
