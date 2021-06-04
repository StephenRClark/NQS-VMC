// ------------------------------------------------------------------------------------------------------------
// Hamiltonian.cpp - contains the declaration of the Hamiltonian object class, which acts as a container of 
//					 Operators and energy terms.
// ------------------------------------------------------------------------------------------------------------

#ifndef HAMILTONIAN_HH
#define HAMILTONIAN_HH

#include <Eigen/Dense>
#include <vector>
#include "Hilbert/Hilbert.hh"
#include "Ansatz/Ansatz.hh"
#include "Operator/Operator.hh"

namespace nqsvmc
{
	class Hamiltonian
	{
	private:
		vector<double> HParams; // Energy terms in the Hamiltonian.
		vector<Operator*> op_; // Vector of Operator pointers.
	public:
		// Constructor requires passing of HParams vector and vector of pointers.
		Hamiltonian(vector<double> hp, vector<Operator*> opptrs);
		double EnergySample(Config* Cfg, Ansatz* AnsObj);
	};
}

#endif