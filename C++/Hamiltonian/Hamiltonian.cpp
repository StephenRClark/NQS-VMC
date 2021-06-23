// ------------------------------------------------------------------------------------------------------------
// Hamiltonian.cpp - contains the implementation of the Hamiltonian object class, which acts as a container of 
//					 Operators and energy terms.
// ------------------------------------------------------------------------------------------------------------

#ifndef HAMILTONIAN_CPP
#define HAMILTONIAN_CPP

#include <Eigen/Dense>
#include <vector>
#include "Hamiltonian.hh"

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	// Definitions of functions in Hamiltonian:
	// Constructor for the Hamiltonian object from parameters and Operators.
	Hamiltonian::Hamiltonian(vector<double> hp, vector<Operator*> opptrs)
	{
		if (hp.size() != opptrs.size())
		{
			cerr << "Number of energy terms does not match number of Operators." << endl;
			std::abort();
		}
		HParams = hp;
		op_ = opptrs;
	}
	// EnergySample gives the local energy from the provided Ansatz and Configuration.
	double Hamiltonian::EnergySample(Config* Cfg, Ansatz* AnsObj)
	{
		double EnLoc = 0;
		MatrixXd OpSamp;
		for (int h = 0; h < HParams.size(); h++)
		{
			OpSamp = op_[h]->LocalSample(Cfg, AnsObj);
			double OpVal = OpSamp.sum();
			if (isinf(OpVal) || isnan(OpVal))
			{
				OpVal = 0;
			}
			EnLoc += (HParams[h] * OpVal);
		}
		return EnLoc;
	}
}

#endif