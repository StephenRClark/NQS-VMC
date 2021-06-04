// ------------------------------------------------------------------------------------------------------------
// MCSampler.hh - contains the declarations of the MCSampler object class, which executes Markov Chain Monte
//				  Carlo sampling with a given Ansatz and Hamiltonian.
// ------------------------------------------------------------------------------------------------------------

#ifndef MCSAMPLER_HH
#define MCSAMPLER_HH

#include <Eigen/Dense>
#include <vector>
#include <random>
#include "Ansatz/Ansatz.hh"
#include "Hilbert/Hilbert.hh"
#include "Operator/Operator.hh"
#include "Hamiltonian/Hamiltonian.hh"

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	// Forward declaration of the MCSampler class.
	class MCSampler
	{
	private:
		Hilbert* hb_; // Pointer to a Hilbert object for configuration generation / manipulation.
		Hamiltonian* hm_; // Pointer to a Hamiltonian object for energy calculation.
		int Nsamp = 5000; // Number of samples to take.
		int Nequil = 5000; // Number of burn in equilibration steps.
		int Nblock = 10; // Number of Markov chain steps between samples. 
	public:
		// Constructor functions:
		// - Basic usage, only assign necessary fields and use defaults for the rest.
		MCSampler(Hilbert* hlb, Hamiltonian* hmt);
		// - Advanced usage, assign all parameters at initialisation.
		MCSampler(Hilbert* hlb, Hamiltonian* hmt, int Nsample, int Nequilibrate, int Nbl);
		// Observer functions:
		// - NumSamp returns the number of sample steps.
		int NumSamp() const;
		// - NumBlock returns the number of steps between samples.
		int NumBlock() const;

		// Parameter assignment functions:
		// - SetNsamp is used to modify the number of samples.
		void SetNsamp(int Nnew);
		// - SetNblock is used to modify the number of steps per sample.
		void SetNblock(int Nnew);
		// - SetNequil is used to modify the number of equilibration steps before sampling.
		void SetNequil(int Nnew);
		// - SetHilbert will replace the Hilbert pointer with a new one.
		void SetHilbert(Hilbert* hnew);
		// - SetHamiltonian will replace the Hamiltonian pointer with a new one.
		void SetHamiltonian(Hamiltonian* hnew);
		// Sampling functions:
		// - MCMCSample - sample the Ansatz object and output local energy, derivatives and move rate.
		void MCMCSample(Ansatz* AnsObj, double& EnAvg, VectorXd& dP, vector<Operator*> opptrs_,
			vector<MatrixXd>& vals, double& MRate) const;
		// - EvalSample - sample the Ansatz object and output local energy.
		void EvalSample(Ansatz* AnsObj, double& EnAvg, vector<Operator*> opptrs_, vector<MatrixXd>& vals) const;
	};
}

#endif