// ------------------------------------------------------------------------------------------------------------
// OperatorDg.hh - contains the definitions of the OperatorDg subclass, which includes operators that are 
//				   diagonal in the configuration basis.
// ------------------------------------------------------------------------------------------------------------

#ifndef OPERATOR_2SITE_HH
#define OPERATOR_2SITE_HH

#include "Operator.hh"

namespace nqsvmc
{
	// Definition of virtual two-site Operator subclass:
	class GeneralOperator2S : public GeneralOperator
	{
	public:
		virtual MatrixXd LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const = 0;
		virtual MatrixXd LocalSample(Config* Cfg, Ansatz* AnsObj) const = 0;
		virtual MatrixXd GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const = 0;
	};
}

#include "BoseOp2S.hh"

#endif