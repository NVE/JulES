## JulES

JulES is a energy market simulation model for long term planning that uses [TuLiPa](https://github.com/NVE/TuLiPa/) as building blocks. JulES clears the power market in small steps with some of the bids generated from deterministic price prognosis models, and stochastic subsystem models for storage valuation. JulES is designed for fast parallel calculations of water values in each reservoir (or battery storage values), without the need for advanced end values.

This is a prototype to test out ideas and inspire the field. Feedback, ideas and contributions are welcome.

### Motivation:
With the transition towards a renewable-based power system, we need models that can represent the new technologies, markets and dynamics. One of [NVE](https://www.nve.no/english/)’s initiatives to improve our understanding of power market modelling, is the ongoing research project called “Power market modelling in Julia”. The goal of this project is to make an algorithm for simulating the Northern European power market with high temporal resolution, detailed hydropower, uncertainty in weather, and using only open-source software. We want to find out if breaking up the problem into many smaller LP-problems, solving many of them deterministically and with open-source solvers, can give fast and good results.

### JulES

#### Simulation concept

JulES clears the power market in small steps with some of the bids generated from deterministic price prognosis models, and stochastic subsystem models for storage valuation. Each subsystem (e.g. battery or watercourse) can have its own storage valuation model that is adapted to the storage system they represent, and can calculate storage values for each individual storage without the need for advanced end values. This should scale well when we want to solve a complex and detailed power market problem since the smaller subsystems can be solved fast and in parallel. 

The simulation model uses a rolling horizon approach where the underlying models are solved for each time step:
1.	Price prognosis models with different degrees of details (e.g. long-term deterministic aggregated power market model)
2.	Models for valuating storage capacity (e.g. medium- or short-term stochastic hydropower scheduling of individual watercourses based on prices from 1)
3.	Market clearing model (e.g. deterministic bid optimization, with part of the bids from 2)

#### TuLiPa
TODO: Benefits of [TuLiPa](https://github.com/NVE/TuLiPa/)

#### Get an overview of JulES:
- src/prognosis.jl – Code for price prognosis models
- src/stochastic.jl – Code for stochastic sub system models
- src/clearing.jl - Code for market clearing problem
- src/util.jl - Various useful functions

#### See also demos:
- demos/Prototype JulES v1

#### Possible improvements to JulES:
See file "Possible improvements to JulES"

### Contact:
Julien Cabrol: jgrc@nve.no

Harald Endresen:

### Licensing:
Copyright 2023 The Norwegian Water Resources and Energy Directorate, and contributors.

Specific license not decided, but the code will be open-source.
