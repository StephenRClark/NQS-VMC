// ------------------------------------------------------------------------------------------------------------
// OperatorDg.cpp - contains the functions of the two siteOperator subclass. Operators provide the methods for
//					calculating physical properties when provided a configuration and a wavefunction.
// ------------------------------------------------------------------------------------------------------------

#ifndef OPERATOR2S_CPP
#define OPERATOR2S_CPP

#include "Operator2S.hh"

namespace nqsvmc
{
	// <<<<<<------------ BiBj functions ------------>>>>>>
	// Constructor only needs to pass a Graph pointer and a Hilbert pointer.
	BpBm_Bose_Op2S::BpBm_Bose_Op2S(Hilbert* hlb, Graph* grp)
	{
		g_ = grp;
		h_ = hlb;
	}
	// Evaluate the operator for only primary lookup lists and sum the contributions.
	MatrixXd BpBm_Bose_Op2S::LocalSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
	{
		vector<int> cfg_vec = Cfg->FullCfg();
		int N = Cfg->Size();
		int Nbond = g_->Nvecs();
		vector<vector<int>> Bonds = g_->BondRead();
		Diff diff; // Generate Diffs and evaluate in sequence rather than creating a struct array.
		MatrixXd value = MatrixXd::Zero(1, 1);
		double LocVal = 0;
		for (int b = 0; b < Nbond; b++)
		{
			LocVal += Evaluate(cfg_vec, Bonds[b], AnsObj);
		}
		value(0) = LocVal;
		return value;
	}
	MatrixXd BpBm_Bose_Op2S::LocalSample(Config* Cfg, Ansatz* AnsObj) const
	{
		vector<int> cfg_vec = Cfg->FullCfg();
		int N = Cfg->Size();
		int Nbond = g_->Nvecs();
		vector<vector<int>> Bonds = g_->BondRead();
		Diff diff; // Generate Diffs and evaluate in sequence rather than creating a struct array.
		MatrixXd value = MatrixXd::Zero(Nbond, 1);
		double LocVal = 0;
		for (int b = 0; b < Nbond; b++)
		{
			value(b, 0) = Evaluate(cfg_vec, Bonds[b], AnsObj);
		}
		return value;
	}
	// Evaluate the operator for all available lookup lists and sum the contributions.
	MatrixXd BpBm_Bose_Op2S::GraphSample(Config* Cfg, double EnLoc, VectorXd dLogp, Ansatz* AnsObj) const
	{
		vector<int> cfg_vec = Cfg->FullCfg();
		int Ntr = g_->Ntranslate();
		MatrixXd value;
		value.resize(Ntr, 1);
		for (int b = 0; b < Ntr; b++)
		{
			vector<int> Bond = g_->BondSearch(b);
			value(b, 0) += Evaluate(cfg_vec, Bond, AnsObj);
		}
		return value;
	}
	double BpBm_Bose_Op2S::Evaluate(vector<int> cfg_vec, vector<int> list, Ansatz* AnsObj) const
	{
		int N = (int)cfg_vec.size();
		vector<double> cfg_double(cfg_vec.begin(), cfg_vec.end());
		double value = 0;
		double locval = 0;
		double Ratio = 1;
		int nmax = h_->SiteDim() - 1;
		int m;
		Diff diff;
		diff.sign = 1;
		for (int n = 0; n < N; n++)
		{
			m = list[n];
			if (m == n) // Operator becomes number counter.
			{
				value += cfg_double[n];
			}
			else
			{
				diff.num = 2;
				diff.val.resize(2);
				diff.pos.resize(2);
				// First move n -> m. Second bracket checks if new configuration is permitted.
				locval = sqrt((cfg_double[n] * (cfg_double[m] + 1))) * (1 - floor(cfg_double[m] / nmax));
				if (locval > 0) // Local value non-negative, will skip if zero.
				{
					diff.val[0] = -1;
					diff.pos[0] = n;
					diff.val[1] = 1;
					diff.pos[1] = m;
					Ratio = AnsObj->PsiRatio(diff);
					value += (locval * Ratio);
				}
				// Conjugate move m -> n.
				locval = sqrt((cfg_double[m] * (cfg_double[n] + 1))) * (1 - floor(cfg_double[n] / nmax));
				if (locval > 0) // Local value non-negative, will skip if zero.
				{
					diff.val[0] = -1;
					diff.pos[0] = m;
					diff.val[1] = 1;
					diff.pos[1] = n;
					Ratio = AnsObj->PsiRatio(diff);
					value += (locval * Ratio);
				}
			}
		}
		// Forward error correcting - NaN and Inf results removed.
		if (isnan(value) || isinf(value))
		{
			value = 0;
		}
		return value;
	}
}

#endif