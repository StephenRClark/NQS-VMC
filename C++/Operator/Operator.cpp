// ------------------------------------------------------------------------------------------------------------
// Operator.cpp - contains the functions of the Operator object class. Operators provide the methods for
//				  calculating physical properties when provided a configuration and a wavefunction.
// ------------------------------------------------------------------------------------------------------------

#ifndef OPERATOR_CPP
#define OPERATOR_CPP

#include <iostream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "Operator.hh"

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	Operator::Operator(Hilbert* hlb, Graph* grp, const char* handle)
	{
		if (strcmp(handle, "EdP") == 0) // Handle for energy gradient operator.
		{
			o_ = new EneLogDeriv_OpDg();
		}
		else if (strcmp(handle, "dPdQ") == 0)
		{
			o_ = new LogLogDeriv_OpDg();
		}
		else if (strcmp(handle, "Ensq") == 0)
		{
			o_ = new EneSq_OpDg();
		}
		else // Will check handle options according to Hilbert typeID:
		{
			if (hlb->Type() == 'b') // Boson Hilbert.
			{
				if (strcmp(handle, "BpBm") == 0)
				{
					o_ = new BpBm_Bose_Op2S(hlb, grp);
				}
				else if (strcmp(handle, "NiNj") == 0)
				{
					o_ = new NiNj_Bose_OpDg(hlb, grp);
				}
				else if (strcmp(handle, "VarN") == 0)
				{
					o_ = new VarN_Bose_OpDg(hlb, grp);
				}
				else if (strcmp(handle, "DbHl") == 0)
				{
					o_ = new DbHl_Bose_OpDg(hlb, grp);
				}
				else if (strcmp(handle, "OcFr") == 0)
				{
					o_ = new OcFr_Bose_OpDg(hlb, grp);
				}
				else
				{
					cerr << "Invalid operator handle - available operators for Bose Hilbert are 'BpBm' and 'NiNj'." << endl;
					std::abort();
				}
			}
		}
	}
	Operator::Operator(const char* handle) // Second version for energy gradient operators.
	{
		if (strcmp(handle, "EdP") == 0) // Handle for energy gradient operator.
		{
			o_ = new EneLogDeriv_OpDg();
		}
		else if (strcmp(handle, "dPdQ") == 0)
		{
			o_ = new LogLogDeriv_OpDg();
		}
		else
		{
			cerr << "Handle-only construction only permitted for energy gradient operators 'EdP' and 'dPdQ'." << endl;
			std::abort();
		}
	}
	// LocalSample will output the local value of the operator evaluated on a particular configuration.
	// - Regardless of single value or multiple-value, output will be a vector for simplicity.
	MatrixXd Operator::LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
	{
		return o_->LocalSample(Cfg, EnLoc, dLogp, AnsObj);
	}
	// - Two versions - one with calls of EnLoc and dLogp, one without.
	MatrixXd Operator::LocalSample(Config* Cfg, Ansatz* AnsObj) const
	{
		return o_->LocalSample(Cfg, AnsObj);
	}
	// GraphSample will output values of the operator for all lookup lists in the associated Graph.
	MatrixXd Operator::GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
	{
		return o_->GraphSample(Cfg, EnLoc, dLogp, AnsObj);
	}
}
#endif