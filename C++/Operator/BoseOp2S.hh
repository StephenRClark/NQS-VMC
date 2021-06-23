// ------------------------------------------------------------------------------------------------------------
// BoseOp2S.hh - contains the definitions of Operator2S objects intended for use with bosonic systems.
// ------------------------------------------------------------------------------------------------------------

#ifndef BOSE_OP_2S_HH
#define BOSE_OP_2S_HH

#include <Eigen/Dense>
#include <vector>
#include "Ansatz/Ansatz.hh"
#include "Graph/Graph.hh"
#include "Hilbert/Hilbert.hh"
#include "Operator2S.hh"

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	class BpBm_Bose_Op2S : public GeneralOperator2S
	{
	private:
		Graph* g_; // Pointer to a Graph for lookup tables. 
		Hilbert* h_; // Pointer to a Hilbert for maximum occupation limits.
	public:
		// Constructor only needs to pass a Graph pointer and a Hilbert pointer.
		BpBm_Bose_Op2S(Hilbert* hlb, Graph* grp);
		// Evaluate the operator for only primary lookup lists and sum the contributions.
		MatrixXd LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const;
		MatrixXd LocalSample(Config* Cfg, Ansatz* AnsObj) const;
		// Evaluate the operator for all available lookup lists and sum the contributions.
		MatrixXd GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const;
		double Evaluate(vector<int> cfg_vec, vector<int> list, Ansatz* AnsObj) const;
	};
}

#endif