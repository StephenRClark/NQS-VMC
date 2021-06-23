// ------------------------------------------------------------------------------------------------------------
// EnGradOpDg.hh - contains the definitions of OperatorDg objects intended for calculating energy gradients.
// ------------------------------------------------------------------------------------------------------------

#ifndef ENEGRAD_OP_DIAG_HH
#define ENEGRAD_OP_DIAG_HH

#include <Eigen/Dense>
#include <vector>
#include "Ansatz/Ansatz.hh"
#include "Graph/Graph.hh"
#include "Hilbert/Hilbert.hh"
#include "OperatorDg.hh"

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	class EneLogDeriv_OpDg : public GeneralOperatorDg // Used for stochastic reconfiguration.
	{
	public:
		// No properties - only need an empty constructor.
		EneLogDeriv_OpDg();
		MatrixXd LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const;
		MatrixXd LocalSample(Config* Cfg, Ansatz* AnsObj) const; // This probably won't ever get called.
		// Shouldn't need to invoke GraphSample, but in the event it is, it returns the same result.
		MatrixXd GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const;
	};

	class LogLogDeriv_OpDg : public GeneralOperatorDg // Used for stochastic reconfiguration.
	{
	public:
		// No properties - only need an empty constructor.
		LogLogDeriv_OpDg();
		MatrixXd LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const;
		MatrixXd LocalSample(Config* Cfg, Ansatz* AnsObj) const; // This probably won't ever get called.
		// Shouldn't need to invoke GraphSample, but in the event it is, it returns the same result.
		MatrixXd GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const;
	};

	class EneSq_OpDg : public GeneralOperatorDg // Used for stochastic reconfiguration.
	{
	public:
		// No properties - only need an empty constructor.
		EneSq_OpDg();
		MatrixXd LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const;
		MatrixXd LocalSample(Config* Cfg, Ansatz* AnsObj) const; // This probably won't ever get called.
		// Primary use case - only used for EvalSample.
		MatrixXd GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const;
	};
}

#endif