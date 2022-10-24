# NQS-VMC
Neural-network Quantum State - Variational Quantum Monte Carlo: a MATLAB implementation
===================================================================================================

Copyright (C) 2019 **Michael Pei** and **Stephen Clark**

![Image](../../images/nqs_vmc.png "nqs_vmc")

Modifier
--------------------

This page documents the details of the Modifier class and related functions. 

N.B: This documentation is a work in progress.

[Back to index.](MATLAB/Documentation/index)

Classes
--------------------

* `Modifier`
    * Abstract class, can only instantiate with subclasses.      
    * Class properties include:
        * `VFlag`
        * `OptInds`
        * `Np`
        * `Graph`
        * `FullCfg`

* `Gutz` - Modifier subclass

* `Jast` - Modifier subclass

* `NNMB` - Modifier subclass

* `NQS` - Modifier subclass

* `NQSA` - Modifier subclass

* `NQSB` - Modifier subclass

* `NQSM` - Modifier subclass

* `NQSU` - Modifier subclass

* `NQSS1` - Modifier subclass

Functions
--------------------

* `Modifier.RndBatchSelect`
* `Modifier.ChangeFullCfg`
* `Modifier.LogDeriv`
* `Modifier.ParamList`
* `Modifier.ParamLoad`
* `Modifier.PsiCfgUpdate`
* `Modifier.PsiGenerate`
* `Modifier.PsiRatio`
* `Modifier.PsiUpdate`
* `NQS.AddHidden`