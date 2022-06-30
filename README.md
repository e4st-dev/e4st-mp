# E4ST

![GitHub contributors](https://img.shields.io/github/contributors/e4st-dev/e4st-mp?logo=GitHub)
![GitHub last commit](https://img.shields.io/github/last-commit/e4st-dev/e4st-mp/main?logo=GitHub)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)

## Overview
<!-- From Manual -->
The Engineering, Economic, and Environmental Electricity Simulation Tool, or E4ST (pronounced "east"), was developed by faculty and research staff at Cornell, Rensselaer Polytechnic Institute, Arizona State University and at Resources for the Future, with support from the U.S. Department of Energy's CERTS program as well as the Power Systems Engineering Research Center (PSERC). The code itself was developed primarily by Biao Mao, Carlos Murillo-Sanchez and Ray Zimmerman with contributions by others, but the E4ST project as a whole was a much larger enterprise with major contributions and efforts by a larger group.

E4ST is built on top of [MATPOWER](https://matpower.org/), a MATLAB® package for solving power flow and optimal power flow problems.  It consists of a set of software toolboxes that can be used to estimate present and future operating and investment states of an electric power system, including:
* generator dispatches
* generator entry and retirement
* locational prices
* fixed and fuel costs
* air emissions
* environmental damages

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
* Make sure that MATPOWER has been added to your MATLAB path via `addpath(genpath(<path_to_matpower>))`
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
* Make sure that MATPOWER has been added to your MATLAB path via `addpath(genpath(<path_to_matpower>))`
* Add the E4ST lib directory recursively to your MATLAB path via `addpath(genpath('lib'));`
* Choose one of the following test cases to run:
    * `test_ces_input.mat`
    * `test_co2cap_input.mat`
    * `test_co2ecap_input.mat`
    * `test_RPSQ_input.mat`
    <!-- TODO: add description of each of these input files -->
* Load the test case, which contains four variables - `esc`, `mpc`, `opt`, and `result`.  Load the test case with the following line of code:
    * ```load(fullfile(e4st_root, 'Input', 'test_ces_input.mat'));```
* Run E4ST with the `runE4ST` command:
    * ```result = runE4ST(mpc, esc, result, opt);```

## Contact Information
Before using E4ST for any major analysis, please contact Ethan Russell at `erussell@rff.org` for additional information on the model, inputs, and outputs.

<!-- TODO: Add Inputs section to describe esc, mpc, opt, and result -->




