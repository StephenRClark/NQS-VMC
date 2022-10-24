# NQS-VMC
Neural-network Quantum State - Variational Quantum Monte Carlo: a MATLAB implementation
===================================================================================================

Copyright (C) 2019 **Michael Pei** and **Stephen Clark**

![Image](../../images/nqs_vmc.png "nqs_vmc")

Optimiser
--------------------

This page documents the details of the Optimiser class and related functions. 

N.B: This documentation is a work in progress.

[Back to index.](MATLAB/Documentation/index)

Classes
--------------------

* `Optimiser`
    * Abstract class, can only instantiate with subclasses.      
    * Class properties include:
        * `Npass`
        * `Ncore`
        * `ExtraSamp`
        * `dERTol`
        * `dEVTol`
        * `dPTol`
        * `ParamTol`
        * `OFrac`

* `StochasticReconfig`

Functions
--------------------

* `Optimiser.SetEnergyTolerances`
* `Optimiser.SetParamTolerances`
* `Optimiser.SetExtraSamples`
* `Optimiser.SetBatchFraction`