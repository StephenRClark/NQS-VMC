// ------------------------------------------------------------------------------------------------------------
// Optimiser.hh - contains the definitions of the Optimiser object class, which modifies variational parameters
//				  of the provided Ansatz using quantities calculated by Monte Carlo sampling.
// ------------------------------------------------------------------------------------------------------------

#ifndef OPTIMISER_HH
#define OPTIMISER_HH

#include <Eigen/Dense>
#include "Ansatz/Ansatz.cpp"
#include "Sampler/MCSampler.cpp"
#include <vector>
#include <random>

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	class StochReconfig;
}

#include "Optimiser/StochReconfig.hh"

#endif