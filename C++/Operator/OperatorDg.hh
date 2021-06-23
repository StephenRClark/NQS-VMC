// ------------------------------------------------------------------------------------------------------------
// OperatorDg.hh - contains the definitions of the OperatorDg subclass, which includes operators that are 
//				   diagonal in the configuration basis.
// ------------------------------------------------------------------------------------------------------------

#ifndef OPERATOR_DIAG_HH
#define OPERATOR_DIAG_HH

#include "Operator.hh"

namespace nqsvmc
{
	// Definition of virtual diagonal Operator subclass:
	class GeneralOperatorDg : public GeneralOperator
	{
	public:
		virtual MatrixXd LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const = 0;
		virtual MatrixXd LocalSample(Config* Cfg, Ansatz* AnsObj) const = 0;
		virtual MatrixXd GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const = 0;
	};
}

#include "BoseOpDg.hh"
#include "EnGradOpDg.hh"

#endif