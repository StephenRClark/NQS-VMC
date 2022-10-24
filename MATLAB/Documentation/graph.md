# NQS-VMC
Neural-network Quantum State - Variational Quantum Monte Carlo: a MATLAB implementation
===================================================================================================

Copyright (C) 2019 **Michael Pei** and **Stephen Clark**

![Image](../../images/nqs_vmc.png "nqs_vmc")

Graph
--------------------

This page documents the details of the Graph class and related functions. 

N.B: This documentation is a work in progress.

[Back to index.](MATLAB/Documentation/index)

Classes
--------------------

* `Graph` - the object class that contains details of the lattice of the chosen wavefunction / model. This class is used in the initialisation of objects in the Operator, Modifier and Reference classes.
    * Constructed as `Graph(Dim,Bonds,Bound,LVecs,SFlag)`
        * `Dim` - `1 x d` vector, `d` is number of independent dimensions.
        * `Bonds` - `N x b` array, `b` is number of connections per site.
        * `Bound` - `1 x d` vector.
        * `LVecs` - `v x d` array, `v` is number of desired primitive lattice vectors.
        * `SFlag` - 1 allows `BondMap` to be generated, 0 disables this.
    * Class properties include:
        * `N` - number of sites within the lattice. Set at initialisation.
        * `Dim` - independent dimensions of the lattice. For example, `Dim = [4 6]` is a lattice with 4 sites along `x` and 6 sites along `y`.
        * `Bound` - boundary conditions along each specified dimension in `Dim`. 0 corresponds to an open / impassable boundary, 1 corresponds to a periodic boundary.
        * `LVecs` - primitive lattice vectors. Can be used to specify desired translational symmetries.
        * `Bonds` - lists of directly bonded sites within the lattice. `Bonds(n,:)` lists all the sites directly linked to site `n`.
        * `SLInds` - `N x 1` list of sublattice indices by site. Can be programatically populated with `FindSublattice`, is set to ones otherwise.
        * `ExtraLabels` - cell array of additional labels for sites. Empty by default, can be manually assigned.
        * `BondMap` - a cell list of `N x 1` bonded site lists, including indirect bonds through combinations of lattice vectors.
        * `Ntr` - number of unique translations in the lattice.
        * `VecInds` - a list of vector indices corresponding to lists in the BondMap. For example, `VecInds(3,:) = [1 2]` means `BondMap{3}` corresponds to translation by the first vector of `LVecs` and twice the second vector, `T = 1 x LVecs(1,:) + 2 x LVecs(2,:)`.

* `HypCub` - a Graph subclass that describes hypercubic lattices. All properties in common with `Graph`.

Functions
--------------------

* `GraphMap` - generates `BondMap`.
* `FindSublattice` - programmatically calculates `SLInds`.
* `ReverseBond` - creates full list of neighbours by doing reverse searches through `Bonds`.
* `Graph2Array` - creates connectivity matrix from `Bonds`.
* `RotateBonds` - for a 2D Graph, calculates `Bonds` after specified number of 90 degree rotations.
* `EdgeLabel` - for a square dual graph, labels sites by their edge number.
* `SquarePlaqLabel` - for a square dual graph, labels sites by their plaquette number.