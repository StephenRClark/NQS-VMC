// ------------------------------------------------------------------------------------------------------------
// Operator.hh - contains the definitions of the Operator object class. Operators provide the methods for
//				 calculating physical properties when provided a configuration and a wavefunction.
// ------------------------------------------------------------------------------------------------------------

#ifndef OPERATOR_HH
#define OPERATOR_HH

#include <Eigen/Dense>
#include <vector>
#include "Ansatz/Ansatz.hh"
#include "Graph/Graph.hh"
#include "Hilbert/Hilbert.hh"

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	class GeneralOperator // Several operator subclasses differing in methods, will need an implementer subclass.
	{
	public:
		// LocalSample will output the local value of the operator evaluated on a particular configuration.
		// - Regardless of single value or multiple-value, output will be a vector for simplicity.
		virtual MatrixXd LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const = 0;
		// - Two versions - one with calls of EnLoc and dLogp, one without.
		virtual MatrixXd LocalSample(Config* Cfg, Ansatz* AnsObj) const = 0;
		// GraphSample will output values of the operator for all lookup lists in the associated Graph.
		virtual MatrixXd GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const = 0;
		// No implementation of CorrMatEls or RWMatEls - those are only required for Operator composition.
	};

	// Forward declaration of implementer subclass:
	class Operator : public GeneralOperator
	{
	private:
		GeneralOperator* o_;
	public:
		Operator(Hilbert* hlb, Graph* grp, const char* handle);
		Operator(const char* handle); // Second version for energy gradient operators.
		// LocalSample will output the local value of the operator evaluated on a particular configuration.
		// - Regardless of single value or multiple-value, output will be a vector for simplicity.
		MatrixXd LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const;
		// - Two versions - one with calls of EnLoc and dLogp, one without.
		MatrixXd LocalSample(Config* Cfg, Ansatz* AnsObj) const;
		// GraphSample will output values of the operator for all lookup lists in the associated Graph.
		MatrixXd GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const;
	};
}

#include "OperatorDg.hh"
#include "Operator2S.hh"

#endif