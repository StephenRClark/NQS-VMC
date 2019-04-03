# NQS-VMC
Neural-network Quantum State - Variational Quantum Monte Carlo: a MATLAB implementation
===================================================================================================

Copyright (C) 2019 **Michael Pei** and **Stephen Clark**

![Image](./images/nqs_vmc.png "nqs_vmc")

Purpose of code
--------------------
Neural-network Quantum States (NQS) are a powerful ansatz for many-body quantum systems inspired by artificial neural networks approaches in machine learning. They were first introduced into the context of variational Monte Carlo (VMC) calculation by the seminal work of Carleo and Troyer _“Solving the quantum many-body problem with artificial neural networks”_ Science **355**, 602 (2017). **NQS-VMC** is a open-source software project with two key aims:

  1. to provide code for performing cutting-edge VMC calculations with NQS,
  2. to be as simple and easier to use as possible to enhance fundamental and pedagogical understanding of NQS.

We have chosen object-orientated MATLAB code to achieve these aims. 

Running the code
--------------------

Running NQS-VMC requires MATLAB. Download the code from this repository. Edit the file `set_NQS_path.m` by defining the directory where the NQS-VMC is located. Run this script to add this directory and all subdirectories to the MATLAB path.

Start from the files `examples/*.m` where scripts are provided for doing NQS VMC calculations on the Bose and Fermi Hubbard model.

Coming soon ...
--------------------
PDF documentation for the code and a discussion of simple examples is in preparation.
  
Key publications
--------------------
The following papers are related to this project:

  * Stephen R. Clark, _"Unifying Neural-network Quantum States and Correlator Product States via Tensor Networks"_,  J. Phys. A: Math. Theor. **51** 135301 (2018), see also [arXiv preprint](https://arxiv.org/abs/1710.03545).
  * Sarah Al-Assam, Stephen R. Clark, Dieter Jaksch, _"The Tensor Network Theory Library"_, J. Stat. Mech. 093102 (2017), see also [arXiv preprint](https://arxiv.org/abs/1610.02244).
  * Sarah Al-Assam, Stephen R. Clark, Chris J. Foot, Dieter Jaksch, _"Capturing long range correlations in two-dimensional quantum lattice systems using correlator product states"_, Phys. Rev. B **84**, 205108 (2011), see also [arXiv preprint](https://arxiv.org/abs/1107.0936).

Funding
--------------------
This code is a part of the project _"Emerging correlations from strong driving: a tensor network projection variational Monte Carlo approach to 2D quantum lattice systems"_ being undertaken via the PI Dr Stephen Clark at the **University of Bristol**, UK and funded by the UK's Engineering and Physical Sciences Research Council (EPSRC) grant EP/P025110/1 - see [GOW page](https://gow.epsrc.ukri.org/NGBOViewGrant.aspx?GrantRef=EP/P025110/1).

![Image](./images/epsrc.png "epsrc")

License
--------------------
This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


