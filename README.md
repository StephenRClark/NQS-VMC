# NQS-VMC
Neural-network Quantum State - Variational Quantum Monte Carlo: a MATLAB implementation
===================================================================================================

Copyright (C) 2019 Michael Pei and Stephen Clark

Neural-network Quantum States (NQS) are a powerful ansatz for many-body quantum systems inspired by artificial neural networks approaches in machine learning. They were first introduced into the context of variational Monte Carlo (VMC) calculation by the seminal work of Carleo and Troyer _“Solving the quantum many-body problem with artificial neural networks”_ Science **355**, 602 (2017). **NQS-VMC** is a open-source software project with a two-fold aim:

  1. providing code for performing cutting-edge VMC calculations with NQS,
  2. being as simple and easier to use as possible to enhance pedagogical and fundamental understanding of NQS.

We have chosen object-orientated MATLAB code to achieve this aims. 


Funding
--------------------
This code is a part of the project _"Emerging correlations from strong driving: a tensor network projection variational Monte Carlo approach to 2D quantum lattice systems"_ being undertaken at the University of Bristol, UK and funded by the UK's Engineering and Physical Sciences Research Council (EPSRC) grant EP/P025110/1 - see [GOW page](https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/P025110/1).





This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

Running the code
--------------------

NQS-VMC requires Matlab

Download code and 
Start from the files examples/*.py.

    > cd lib
    > mkdir build
    > cd build
    > CXX=icpc CC=icc cmake .. -DMKL=ON
    > make
    > cd ../..


Bullet list:

  * apples
  * oranges
  * pears

Numbered list:

  1. wash
  2. rinse
  3. repeat

A [example].

  [example]: http://example.com

![Image](Icon-pictures.png "icon")



Performance testing
-------------------

### 1. Find the most costly functions:


