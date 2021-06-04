// ------------------------------------------------------------------------------------------------------------
// OperatorDg.hh - contains the definitions of the OperatorDg subclass, which includes operators that are 
//				   diagonal in the configuration basis.
// ------------------------------------------------------------------------------------------------------------

#ifndef OPERATOR_2SITE_HH
#define OPERATOR_2SITE_HH

#include <Eigen/Dense>
#include <vector>
#include "Ansatz/Ansatz.hh"
#include "Graph/Graph.hh"
#include "Hilbert/Hilbert.hh"
#include "Operator.hh"

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	class GeneralOperator2S : public GeneralOperator
	{
	public:
		virtual MatrixXd LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const = 0;
		virtual MatrixXd LocalSample(Config* Cfg, Ansatz* AnsObj) const = 0;
		virtual MatrixXd GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const = 0;
	};

	// Forward declarations of subclasses.
}

#include "BoseOp2S.hh"

#endif