# E4ST

![GitHub contributors](https://img.shields.io/github/contributors/e4st-dev/e4st-mp?logo=GitHub)
![GitHub last commit](https://img.shields.io/github/last-commit/e4st-dev/e4st-mp/main?logo=GitHub)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

## Overview
<!-- From Manual -->
The Engineering, Economic, and Environmental Electricity Simulation Tool (E4ST, pronounced “east”) is policy analysis and planning software built to simulate in detail how the power sector will operate and evolve in response to environmental and non-environmental policies and regulations, renewable and non-renewable generation investments, transmission investments, demand changes, and so on. The user specifies the policies, investments, and other inputs of each simulation. E4ST predicts operation, generator investment and retirement, prices, consumer welfare, producer profits, emissions, and health effects of emissions, among other outcomes, in each simulated year. It is therefore well suited for comprehensive benefit-cost analysis. 

As implemented by its developers at Resources for the Future, E4ST has exceptionally high spatial and engineering detail relative to other US national policy analysis models, using a detailed model of the US power grid with over 5,000 nodes (substations). Within E4ST, power flows on the grid are represented in an unusually realistic manner through linear approximations of the non-linear equations that determine actual transmission line flows. E4ST is also a power system planning model. As a planning model, its distinguishing strengths include endogenous prediction of generator construction and retirement, and the aforementioned comprehensive benefit-cost analysis capabilities.

At the heart of E4ST is an optimization problem that represents the behavior of the system operators, electricity end-users, generators, and generation developers. It is a cost minimization problem that minimizes the sum of generator variable costs, generator fixed costs, end-user consumer surplus losses, and investment costs. This problem simultaneously mimics the decisions of system operators and the investment and retirement decisions of generation owners. 

E4ST is built on top of [MATPOWER](https://matpower.org/), a MATLAB® package for solving power flow and optimal power flow problems.

E4ST has been developed by researchers at Resources for the Future, Cornell University, Arizona State University, and Rensselaer Polytechnic Institute, with support from the U.S. Department of Energy's CERTS program, the Power Systems Engineering Research Center (PSERC), the Sloan Foundation, Breakthrough Energy, the New York Independent System Operator, and other funders of Resources for the Future.

[E4ST.com](https://e4st.com/) is an older website with additional information about E4ST, that complements the information available at this site.


## E4ST Model Formulation
Find the E4ST Model formulation in [`docs/E4ST_Formulation.pdf`](https://github.com/e4st-dev/e4st-mp/blob/main/docs/E4ST_Formulation.pdf)

## Getting Started

### System Requirements
* MATLAB® version R2016b or later
* MATPOWER 7.0 or later
* Optional: Commercial solver Gurobi v9.1.2 or later

### E4ST Installation
E4ST can be installed by either of the following methods:
* Clone the repository using a command such as:
    * `git clone git@github.com:e4st-dev/e4st-mp.git`
* Download the repository as a zip file

### Testing E4ST
To test that your install of E4ST is working correctly, there is a small suite of tests that can be run with the following procedure:
* Make sure that MATPOWER has been added to your MATLAB path via `run('<path to MATPOWER>/install_matpower.m')`
* Add the E4ST lib directory recursively to your MATLAB path via `addpath(genpath('lib'));`
* Run `test_e4st`
    * This should run without error and print something similar to the following:
        ```
        t_apply_changes....ok
        t_e4st_solve.......ok
        t_e4st_caplim......ok
        t_e4st_storage.....ok
        t_e4st_dac.........ok
        All tests successful (1249 of 1249)
        Elapsed time 6.62 seconds.
        ```

### Running Sample Cases
* Make sure that MATPOWER has been added to your MATLAB path via `run('<path to MATPOWER>/install_matpower.m')`
* Add the E4ST lib directory recursively to your MATLAB path via `addpath(genpath('lib'));`
* Choose one of the following test cases to run. Each test case contains the inputs for a 30-bus power system, with varying power sector policies and generator details:
    * `test_ces_input.mat`
    * `test_co2cap_input.mat`
    * `test_co2ecap_input.mat`
    * `test_RPSQ_input.mat`
    <!-- TODO: add description of each of these input files -->
* Load the test case, which contains four variables - `esc`, `mpc`, `opt`, and `result`.  Load the test case with the following line of code:
    * ```load(fullfile(e4st_root, 'Input', 'test_ces_input.mat'));```
* Run E4ST with the `RunE4ST` command:
    * ```result = RunE4ST(mpc, esc, result, opt);```

## Contact Information
Before using E4ST for any major analysis, please contact Ethan Russell at `erussell@rff.org` for additional information on the model, inputs, and outputs.

<!-- TODO: Add Inputs section to describe esc, mpc, opt, and result -->




