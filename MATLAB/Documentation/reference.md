# NQS-VMC
Neural-network Quantum State - Variational Quantum Monte Carlo: a MATLAB implementation
===================================================================================================

Copyright (C) 2019 **Michael Pei** and **Stephen Clark**

![Image](../../images/nqs_vmc.png "nqs_vmc")

Reference
--------------------

This page documents the details of the Reference class and related functions. 

N.B: This documentation is a work in progress.

[Back to index.](MATLAB/Documentation/index)

Classes
--------------------

* `Reference`
    * Abstract class, can only instantiate with subclasses.      
    * Class properties include:
        * `VFlag`
        * `Type`
        * `Np`
        * `FullCfg`

* `BECR` - Reference subclass
    * Subclass specific properties include:
        * `Nb`
        * `Np`
        * `Graph`
        * `Occ`
        * `SPO`

* `Plus` - Reference subclass

* `Pfaf` - Reference subclass

* `SDet` - Reference subclass

Functions
--------------------

* `Reference.PrepPsi`
* `Reference.PsiCfgUpdate`
* `Reference.PsiGenerate`
* `Reference.PsiRatio`