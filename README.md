# NQS-VMC
Neural-network Quantum State - Variational Quantum Monte Carlo: a Matlab implementation
===================================================================================================

Copyright (C) 2019 Michael Pei and Stephen Clark

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



Performance testing
-------------------

### 1. Find the most costly functions:


