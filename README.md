Sample codes for the paper "Revealing patterns of cultural transmission from frequency data: equilibrium  and non-equilibrium assumptions". Notice that in order to reproduce the results of the paper the simulation needs to be executed on a High Performance Cluster. The examples script provides the basic workflow with just 0.01% of the runs used in the paper.

## Paper Reference 
Crema,E.R, Kandler, A., Shennan, S.J. "Revealing patterns of cultural transmission from frequency data:
equilibrium and non-equilibrium assumptions." Sci. Rep. 6, 39122; doi: 10.1038/srep39122 (2016).

## Author of the Repository:
Enrico R. Crema (enrico.crema@gmail.com)

## Contents
* ./src/ ... contains core simulation model for the equilibrium, variable population, and variable population/transmission versions as well as a general utility function.
* ./data/observedFrequencies.csv ... frequencies of different decorative motifs (coded BT1, BT2, ...) for phases VIII to XIV in the Merzbach assemblage
* ./logEquilibrium.R, ./logVarPopTrans.R, ./logVarPop.R ... sample script for executing the ABC framework. 

## Licences
Text: CC-BY (http://creativecommons.org/licenses/by/4.0/)

Code: MIT (http://opensource.org/licenses/MIT year: 2014, copyright holder: Enrico R. Crema)

Data: CC0 (http://creativecommons.org/publicdomain/zero/1.0/) attribution requested in reuse

## Dependencies
R version 3.3.0 (2016-05-03)
Platform: x86_64-apple-darwin13.4.0 (64-bit)

attached base packages:
* stats
* graphics
* grDevices
* utils
* datasets
* methods
* base    

other attached packages:
* LaplacesDemon_15.03.19

loaded via a namespace (and not attached):
* parallel_3.3.0
