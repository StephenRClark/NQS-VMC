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
		NiNj_Bose_OpDg(Hilbert* hlb, Graph* grp)
		{
			g_ = grp;
		}
		// Evaluate the operator for only primary lookup lists and sum the contributions.
		MatrixXd LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
		{
			int Nbond = g_->Nvecs();
			vector<vector<int>> Bonds = g_->BondRead();
			vector<int> cfg_vec = Cfg->FullCfg();
			MatrixXd value = MatrixXd::Zero(1, 1);
			for (int b = 0; b < Nbond; b++)
			{
				value(0) += Evaluate(cfg_vec, Bonds[b]);
			}
			return value;
		}
		MatrixXd LocalSample(Config* Cfg, Ansatz* AnsObj) const
		{
			int Nbond = g_->Nvecs();
			vector<vector<int>> Bonds = g_->BondRead();
			vector<int> cfg_vec = Cfg->FullCfg();
			MatrixXd value = MatrixXd::Zero(Nbond, 1);
			for (int b = 0; b < Nbond; b++)
			{
				value(b,0) += Evaluate(cfg_vec, Bonds[b]);
			}
			return value;
		}
		// Evaluate the operator for all available lookup lists and sum the contributions.
		MatrixXd GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
		{
			int Ntr = g_->Ntranslate();
			vector<int> cfg_vec = Cfg->FullCfg();
			MatrixXd value;
			value.resize(Ntr, 1);
			for (int n = 0; n < Ntr; n++)
			{
				vector<int> Bond = g_->BondSearch(n);				
				value(n, 0) = Evaluate(cfg_vec, Bond);
			}
			return value;
		}
		double Evaluate(vector<int> cfg_vec, vector<int> list) const
		{
			int N = (int)cfg_vec.size();
			double offset = 0;
			double value = 0;
			int m = 0;
			if (list[0] == 0)
			{
				offset = 1;
			}
			for (int n = 0; n < N; n++)
			{
				m = list[n];
				value += (cfg_vec[n] * (cfg_vec[m] - offset));
			}
			return value;
		}
	};

	class DbHl_Bose_OpDg : public GeneralOperatorDg
	{
	private:
		Graph* g_; // Pointer to a Graph for lookup tables. 
	public:
		// Constructor only needs to pass a Graph pointer.
		DbHl_Bose_OpDg(Hilbert* hlb, Graph* grp)
		{
			g_ = grp;
		}
		// Evaluate the operator for only primary lookup lists and sum the contributions.
		MatrixXd LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
		{
			int Nbond = g_->Nvecs();
			vector<vector<int>> Bonds = g_->BondRead();
			vector<int> cfg_vec = Cfg->FullCfg();
			MatrixXd value = MatrixXd::Zero(1, Bonds.size());
			for (int b = 0; b < Nbond; b++)
			{
				value(b) += Evaluate(cfg_vec, Bonds[b]);
			}
			return value;
		}
		MatrixXd LocalSample(Config* Cfg, Ansatz* AnsObj) const
		{
			int Nbond = g_->Nvecs();
			vector<vector<int>> Bonds = g_->BondRead();
			vector<int> cfg_vec = Cfg->FullCfg();
			MatrixXd value = MatrixXd::Zero(Nbond, 1);
			for (int b = 0; b < Nbond; b++)
			{
				value(b, 0) += Evaluate(cfg_vec, Bonds[b]);
			}
			return value;
		}
		// Evaluate the operator for all available lookup lists and sum the contributions.
		MatrixXd GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
		{
			int Ntr = g_->Ntranslate();
			vector<int> cfg_vec = Cfg->FullCfg();
			MatrixXd value;
			value.resize(Ntr, 1);
			for (int n = 0; n < Ntr; n++)
			{
				vector<int> Bond = g_->BondSearch(n);
				value(n, 0) = Evaluate(cfg_vec, Bond);
			}
			return value;
		}
		double Evaluate(vector<int> cfg_vec, vector<int> list) const
		{
			int N = (int)cfg_vec.size();
			double value = 0;
			double db_i = 0;
			double hl_i = 0;
			double db_j = 0;
			double hl_j = 0;
			int m = 0;
			for (int n = 0; n < N; n++)
			{
				m = list[n];
				db_i = (double)(cfg_vec[n] == 2);
				hl_i = (double)(cfg_vec[n] == 0);
				db_j = (double)(cfg_vec[m] == 2);
				hl_j = (double)(cfg_vec[m] == 0);				
				value += 0.5*(db_i*hl_j + db_j*hl_i)/N;
			}
			return value;
		}
	};

	class VarN_Bose_OpDg : public GeneralOperatorDg
	{
	private:
		Graph* g_; // Pointer to a Graph for lookup tables. 
	public:
		// Constructor only needs to pass a Graph pointer.
		VarN_Bose_OpDg(Hilbert* hlb, Graph* grp)
		{
			g_ = grp;
		}
		// Evaluate the operator for only primary lookup lists and sum the contributions.
		MatrixXd LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
		{
			vector<int> cfg_vec = Cfg->FullCfg();
			int Nb = Cfg->Nparticle();
			MatrixXd value = MatrixXd::Zero(1, 1);
			value(0) = Evaluate(cfg_vec, Nb);
			return value;
		}
		MatrixXd LocalSample(Config* Cfg, Ansatz* AnsObj) const
		{
			vector<int> cfg_vec = Cfg->FullCfg();
			int Nb = Cfg->Nparticle();
			MatrixXd value = MatrixXd::Zero(1, 1);
			value(0) = Evaluate(cfg_vec, Nb);
			return value;
		}
		// Evaluate the operator for all available lookup lists and sum the contributions.
		MatrixXd GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
		{
			vector<int> cfg_vec = Cfg->FullCfg();
			int Nb = Cfg->Nparticle();
			MatrixXd value = MatrixXd::Zero(1, 1);
			value(0) = Evaluate(cfg_vec, Nb);
			return value;
		}
		double Evaluate(vector<int> cfg_vec, int Nb) const
		{
			int N = (int)cfg_vec.size();
			double value = 0;
			double cn = 0;
			for (int n = 0; n < N; n++)
			{
				cn = (double)cfg_vec[n];
				value += (pow(cn - ((double)Nb)/((double)N),2)/((double)N));
			}
			return value;
		}
	};

	class OcFr_Bose_OpDg : public GeneralOperatorDg
	{
	private:
		Graph* g_; // Pointer to a Graph for lookup tables. 
		int dim; // On-site dimension - need to record from initialisation Hilbert.
	public:
		// Constructor only needs to pass a Graph pointer.
		OcFr_Bose_OpDg(Hilbert* hlb, Graph* grp)
		{
			g_ = grp;
			dim = hlb->SiteDim();
		}
		// Evaluate the operator for only primary lookup lists and sum the contributions.
		MatrixXd LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
		{
			vector<int> cfg_vec = Cfg->FullCfg();
			MatrixXd value = MatrixXd::Zero(1, dim);
			for (int d = 0; d < dim; d++)
			{
				value(0, d) = Evaluate(cfg_vec, d);
			}
			return value;
		}
		MatrixXd LocalSample(Config* Cfg, Ansatz* AnsObj) const
		{
			vector<int> cfg_vec = Cfg->FullCfg();
			MatrixXd value = MatrixXd::Zero(1, dim);
			for (int d = 0; d < dim; d++)
			{
				value(0, d) = Evaluate(cfg_vec, d);
			}
			return value;
		}
		// Evaluate the operator for all available lookup lists and sum the contributions.
		MatrixXd GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
		{
			vector<int> cfg_vec = Cfg->FullCfg();
			MatrixXd value = MatrixXd::Zero(1, dim);
			for (int d = 0; d < dim; d++)
			{
				value(0, d) = Evaluate(cfg_vec, d);
			}
			return value;
		}
		double Evaluate(vector<int> cfg_vec, int d) const
		{
			int N = (int)cfg_vec.size();
			double value = 0;
			for (int n = 0; n < N; n++)
			{
				value += ((double)(cfg_vec[n] == d) / ((double)N));
			}
			return value;
		}
	};
}

#endif