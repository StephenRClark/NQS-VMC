# NQS-VMC
Neural-network Quantum State - Variational Quantum Monte Carlo: a MATLAB implementation
===================================================================================================

Copyright (C) 2019 **Michael Pei** and **Stephen Clark**

![Image](../../images/nqs_vmc.png "nqs_vmc")

Sampler
--------------------

This page documents the details of the Sampler class and related functions. 

N.B: This documentation is a work in progress.

[Back to index.](MATLAB/Documentation/index)

Classes
--------------------

* `Sampler`     
    * Class properties include:
        * `Nequil`
        * `Nblock`
        * `Nsamp`
        * `Hamiltonian`
        * `Operators`

Functions
--------------------

* `Sampler.SetNequil`
* `Sampler.SetNsamp`
* `Sampler.SetNblock`
* `Sampler.SetHamiltonian`
* `Sampler.SetOperator`
* `Sampler.EvalSample`
* `Sampler.MCMCSample`
* `Sampler.MultiChainSample`
