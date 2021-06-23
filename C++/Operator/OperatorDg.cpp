// ------------------------------------------------------------------------------------------------------------
// OperatorDg.cpp - contains the functions of the diagonal Operator subclass. Operators provide the methods for
//					calculating physical properties when provided a configuration and a wavefunction.
// ------------------------------------------------------------------------------------------------------------

#ifndef OPERATORDG_CPP
#define OPERATORDG_CPP

#include "OperatorDg.hh"

namespace nqsvmc
{
	// <<<<<<------------ EdP functions ------------>>>>>>
	// No properties - only need an empty constructor.
	EneLogDeriv_OpDg::EneLogDeriv_OpDg()
	{ }
	MatrixXd EneLogDeriv_OpDg::LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
	{
		VectorXd EdP = dLogp * EnLoc;
		return EdP;
	}
	MatrixXd EneLogDeriv_OpDg::LocalSample(Config* Cfg, Ansatz* AnsObj) const // This probably won't ever get called.
	{
		VectorXd EdP = VectorXd::Zero(1);
		return EdP;
	}
	// Shouldn't need to invoke GraphSample, but in the event it is, it returns the same result.
	MatrixXd EneLogDeriv_OpDg::GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
	{
		VectorXd EdP = dLogp * EnLoc;
		return EdP;
	}

	// <<<<<<------------ dPdQ functions ------------>>>>>>
	// No properties - only need an empty constructor.
	LogLogDeriv_OpDg::LogLogDeriv_OpDg()
	{ }
	MatrixXd LogLogDeriv_OpDg::LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
	{
		MatrixXd dPdQ = dLogp * dLogp.transpose();
		return dPdQ;
	}
	MatrixXd LogLogDeriv_OpDg::LocalSample(Config* Cfg, Ansatz* AnsObj) const // This probably won't ever get called.
	{
		VectorXd dPdQ = VectorXd::Zero(1);
		return dPdQ;
	}
	// Shouldn't need to invoke GraphSample, but in the event it is, it returns the same result.
	MatrixXd LogLogDeriv_OpDg::GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
	{
		MatrixXd dPdQ = dLogp * dLogp.transpose();
		return dPdQ;
	}

	// <<<<<<------------ EnSq functions ------------>>>>>>
	// No properties - only need an empty constructor.
	EneSq_OpDg::EneSq_OpDg()
	{ }
	MatrixXd EneSq_OpDg::LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
	{
		MatrixXd Esq = MatrixXd::Zero(1, 1);
		Esq(0) = pow(EnLoc, 2);
		return Esq;
	}
	MatrixXd EneSq_OpDg::LocalSample(Config* Cfg, Ansatz* AnsObj) const // This probably won't ever get called.
	{
		MatrixXd Esq = MatrixXd::Zero(1, 1);
		return Esq;
	}
	// Primary use case - only used for EvalSample.
	MatrixXd EneSq_OpDg::GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
	{
		MatrixXd Esq = MatrixXd::Zero(1, 1);
		Esq(0) = pow(EnLoc, 2);
		return Esq;
	}

	// <<<<<<------------ NiNj Bose functions ------------>>>>>>
	// Constructor only needs to pass a Graph pointer.
	NiNj_Bose_OpDg::NiNj_Bose_OpDg(Hilbert* hlb, Graph* grp)
	{
		g_ = grp;
	}
	// Evaluate the operator for only primary lookup lists and sum the contributions.
	MatrixXd NiNj_Bose_OpDg::LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
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
	MatrixXd NiNj_Bose_OpDg::LocalSample(Config* Cfg, Ansatz* AnsObj) const
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
	MatrixXd NiNj_Bose_OpDg::GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
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
	double NiNj_Bose_OpDg::Evaluate(vector<int> cfg_vec, vector<int> list) const
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

	// <<<<<<------------ DbHl Bose functions ------------>>>>>>
	// Constructor only needs to pass a Graph pointer.
	DbHl_Bose_OpDg::DbHl_Bose_OpDg(Hilbert* hlb, Graph* grp)
	{
		g_ = grp;
	}
	// Evaluate the operator for only primary lookup lists and sum the contributions.
	MatrixXd DbHl_Bose_OpDg::LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
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
	MatrixXd DbHl_Bose_OpDg::LocalSample(Config* Cfg, Ansatz* AnsObj) const
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
	MatrixXd DbHl_Bose_OpDg::GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
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
	double DbHl_Bose_OpDg::Evaluate(vector<int> cfg_vec, vector<int> list) const
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
			value += 0.5 * (db_i * hl_j + db_j * hl_i) / N;
		}
		return value;
	}

	// <<<<<<------------ VarN Bose functions ------------>>>>>>
	// Constructor only needs to pass a Graph pointer.
	VarN_Bose_OpDg::VarN_Bose_OpDg(Hilbert* hlb, Graph* grp)
	{
		g_ = grp;
	}
	// Evaluate the operator for only primary lookup lists and sum the contributions.
	MatrixXd VarN_Bose_OpDg::LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
	{
		vector<int> cfg_vec = Cfg->FullCfg();
		int Nb = Cfg->Nparticle();
		MatrixXd value = MatrixXd::Zero(1, 1);
		value(0) = Evaluate(cfg_vec, Nb);
		return value;
	}
	MatrixXd VarN_Bose_OpDg::LocalSample(Config* Cfg, Ansatz* AnsObj) const
	{
		vector<int> cfg_vec = Cfg->FullCfg();
		int Nb = Cfg->Nparticle();
		MatrixXd value = MatrixXd::Zero(1, 1);
		value(0) = Evaluate(cfg_vec, Nb);
		return value;
	}
	// Evaluate the operator for all available lookup lists and sum the contributions.
	MatrixXd VarN_Bose_OpDg::GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
	{
		vector<int> cfg_vec = Cfg->FullCfg();
		int Nb = Cfg->Nparticle();
		MatrixXd value = MatrixXd::Zero(1, 1);
		value(0) = Evaluate(cfg_vec, Nb);
		return value;
	}
	double VarN_Bose_OpDg::Evaluate(vector<int> cfg_vec, int Nb) const
	{
		int N = (int)cfg_vec.size();
		double value = 0;
		double cn = 0;
		for (int n = 0; n < N; n++)
		{
			cn = (double)cfg_vec[n];
			value += (pow(cn - ((double)Nb) / ((double)N), 2) / ((double)N));
		}
		return value;
	}

	// <<<<<<------------ OcFr Bose functions ------------>>>>>>
	// Constructor only needs to pass a Graph pointer.
	OcFr_Bose_OpDg::OcFr_Bose_OpDg(Hilbert* hlb, Graph* grp)
	{
		g_ = grp;
		dim = hlb->SiteDim();
	}
	// Evaluate the operator for only primary lookup lists and sum the contributions.
	MatrixXd OcFr_Bose_OpDg::LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
	{
		vector<int> cfg_vec = Cfg->FullCfg();
		MatrixXd value = MatrixXd::Zero(1, dim);
		for (int d = 0; d < dim; d++)
		{
			value(0, d) = Evaluate(cfg_vec, d);
		}
		return value;
	}
	MatrixXd OcFr_Bose_OpDg::LocalSample(Config* Cfg, Ansatz* AnsObj) const
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
	MatrixXd OcFr_Bose_OpDg::GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
	{
		vector<int> cfg_vec = Cfg->FullCfg();
		MatrixXd value = MatrixXd::Zero(1, dim);
		for (int d = 0; d < dim; d++)
		{
			value(0, d) = Evaluate(cfg_vec, d);
		}
		return value;
	}
	double OcFr_Bose_OpDg::Evaluate(vector<int> cfg_vec, int d) const
	{
		int N = (int)cfg_vec.size();
		double value = 0;
		for (int n = 0; n < N; n++)
		{
			value += ((double)(cfg_vec[n] == d) / ((double)N));
		}
		return value;
	}
}

#endif