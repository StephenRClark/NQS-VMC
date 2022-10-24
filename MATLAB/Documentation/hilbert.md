# NQS-VMC
Neural-network Quantum State - Variational Quantum Monte Carlo: a MATLAB implementation
===================================================================================================

Copyright (C) 2019 **Michael Pei** and **Stephen Clark**

![Image](../../images/nqs_vmc.png "nqs_vmc")

Hilbert
--------------------

This page documents the details of the Hilbert class and related functions. 

N.B: This documentation is a work in progress.

[Back to index.](MATLAB/Documentation/index)

Classes
--------------------

* `Hilbert`
    * Abstract class, can only instantiate with subclasses.      
    * Class properties include:
        * `Type`
        * `N`
        * `d`
        * `Sector`
        * `PropMoveFunc`
        * `FullCfg`
        * `RandomCfgFunc`
        * `Diff2CfgFunc`
        * `CParams`

* `Bose` - Hilbert subclass

* `Ferm` - Hilbert subclass

* `Spin` - Hilbert subclass

* `SXYZ` - Hilbert subclass

* `Unary` - Hilbert subclass

Functions
--------------------

* `Hilbert.RandomCfg`
* `Hilbert.PropMove`
* `Hilbert.FullCfg`
* `Hilbert.Diff2Cfg`
* `Hilbert.ChangePropMove`
* `Hilbert.ChangeStartCfg`