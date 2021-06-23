// ------------------------------------------------------------------------------------------------------------
// BoseOpDg.hh - contains the definitions of OperatorDg objects intended for use with bosonic systems.
// ------------------------------------------------------------------------------------------------------------

#ifndef BOSE_OP_DIAG_HH
#define BOSE_OP_DIAG_HH

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
	class NiNj_Bose_OpDg : public GeneralOperatorDg
	{
	private:
		Graph* g_; // Pointer to a Graph for lookup tables. 
	public:
		// Constructor only needs to pass a Graph pointer.
		NiNj_Bose_OpDg(Hilbert* hlb, Graph* grp);
		// Evaluate the operator for only primary lookup lists and sum the contributions.
		MatrixXd LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const;
		MatrixXd LocalSample(Config* Cfg, Ansatz* AnsObj) const;
		// Evaluate the operator for all available lookup lists and sum the contributions.
		MatrixXd GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const;
		double Evaluate(vector<int> cfg_vec, vector<int> list) const;
	};

	class DbHl_Bose_OpDg : public GeneralOperatorDg
	{
	private:
		Graph* g_; // Pointer to a Graph for lookup tables. 
	public:
		// Constructor only needs to pass a Graph pointer.
		DbHl_Bose_OpDg(Hilbert* hlb, Graph* grp);
		// Evaluate the operator for only primary lookup lists and sum the contributions.
		MatrixXd LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const;
		MatrixXd LocalSample(Config* Cfg, Ansatz* AnsObj) const;
		// Evaluate the operator for all available lookup lists and sum the contributions.
		MatrixXd GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const;
		double Evaluate(vector<int> cfg_vec, vector<int> list) const;
	};

	class VarN_Bose_OpDg : public GeneralOperatorDg
	{
	private:
		Graph* g_; // Pointer to a Graph for lookup tables. 
	public:
		// Constructor only needs to pass a Graph pointer.
		VarN_Bose_OpDg(Hilbert* hlb, Graph* grp);
		// Evaluate the operator for only primary lookup lists and sum the contributions.
		MatrixXd LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const;
		MatrixXd LocalSample(Config* Cfg, Ansatz* AnsObj) const;
		// Evaluate the operator for all available lookup lists and sum the contributions.
		MatrixXd GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const;
		double Evaluate(vector<int> cfg_vec, int Nb) const;
	};

	class OcFr_Bose_OpDg : public GeneralOperatorDg
	{
	private:
		Graph* g_; // Pointer to a Graph for lookup tables. 
		int dim; // On-site dimension - need to record from initialisation Hilbert.
	public:
		// Constructor only needs to pass a Graph pointer.
		OcFr_Bose_OpDg(Hilbert* hlb, Graph* grp);
		// Evaluate the operator for only primary lookup lists and sum the contributions.
		MatrixXd LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const;
		MatrixXd LocalSample(Config* Cfg, Ansatz* AnsObj) const;
		// Evaluate the operator for all available lookup lists and sum the contributions.
		MatrixXd GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const;
		double Evaluate(vector<int> cfg_vec, int d) const;
	};
}

#endif