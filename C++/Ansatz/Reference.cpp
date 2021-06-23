// ------------------------------------------------------------------------------------------------------------
// Reference.cpp - contains the functions of the Reference object subclasses. References are system specific
//				   wavefunctions to which Modifiers are applied to form a full Ansatz.
// ------------------------------------------------------------------------------------------------------------

#ifndef REFERENCE_CPP
#define REFERENCE_CPP

#include "Reference.hh"

namespace nqsvmc
{
	// <<<<<<------------ BECR functions ------------>>>>>>
	// Constructors for the Bose-Einstein Condensate Reference:
	// - Initialisation with a single particle Hamiltonian matrix.
	BECR::BECR(Hilbert* hlb, MatrixXd SPH)
	{
		h_ = hlb;
		N = hlb->SysSize();
		Nb = N; // Placeholder, initialised later.
		Occ.resize(N);
		OccP.resize(N);
		if (SPH.rows() != SPH.cols())
		{
			cerr << "Please input a real NxN matrix." << endl;
			std::abort();
		}
		SelfAdjointEigenSolver<MatrixXd> eig(SPH);
		SPO = eig.eigenvectors().col(0);
	}
	// - Initialisation with a single particle orbital vector.
	BECR::BECR(Hilbert* hlb, VectorXd Orbital)
	{
		h_ = hlb;
		N = hlb->SysSize();
		Nb = N; // Placeholder, initialised later.
		NbP = N;
		Occ.resize(N);
		OccP.resize(N);
		if (Orbital.size() != N)
		{
			cerr << "Mismatch between Hilbert system size and orbital system size." << endl;
			std::abort();
		}
		SPO = Orbital;
	}
	// Observer functions for extracting information: (Variation disabled - all return zero).
	// - BECR::Nparam will return 0 - variation is disabled.
	int BECR::Nparam() const
	{
		return 0;
	}
	// - BECR::VarFlag will return 0 - variation is disabled.
	bool BECR::VarFlag() const
	{
		return false;
	}
	// - BECR::OptIndRead will return individual parameter flags for variational modification.
	VectorXd BECR::OptIndRead() const
	{
		VectorXd Inds(1);
		Inds(0) = 0;
		return Inds;
	}
	// - BECR::ParamList will return a vector of the variational parameters in the Reference.
	VectorXd BECR::ParamList() const
	{
		VectorXd Params(1);
		Params(0) = 0;
		return Params;
	}

	// Variational modification management functions: (Variation disabled - all return message)
	// - BECR::VarSwitch will alter VFlag if allowed.
	void BECR::VarSwitch()
	{
		cout << "Variation of this Reference is not supported - switching disabled." << endl;
		return;
	}
	// - BECR::OptIndLoad will load a vector of updated parameter flags (values 0/1).
	void BECR::OptIndLoad(VectorXd newinds)
	{
		cout << "Variation of this Reference is not supported - switching disabled." << endl;
		return;
	}
	// - BECR::RndBtcSelect will randomly alter the individual parameter flags.
	// -- Default zero argument version will select 1/e of the parameters.
	void BECR::RndBtcSelect()
	{
		cout << "Variation of this Reference is not supported - switching disabled." << endl;
		return;
	}
	// -- Second version with fraction will select input fraction if valid.
	void BECR::RndBtcSelect(double pfrac)
	{
		cout << "Variation of this Reference is not supported - switching disabled." << endl;
		return;
	}
	// - SetParamCap will change the maximum allowed amplitude of a single parameter.
	void BECR::SetParamCap(double newcap)
	{
		cout << "Variation of this Reference is not supported - parameter capping disabled." << endl;
		return;
	}
	// - ParamLoad will change the existing parameters to the provided parameters.
	void BECR::ParamLoad(VectorXd NewParams)
	{
		cout << "Variation of this Reference is not supported - parameter loading disabled." << endl;
		return;
	}

	// Wavefunction calculation functions:


	// - BECR::CfgRead will allow for additional modification of FullCfg methods in Config.
	VectorXd BECR::CfgRead(Config* Cfg) const
	{
		vector<int> cfg_init = Cfg->FullCfg();
		vector<double> cfg_double(cfg_init.begin(), cfg_init.end());
		double* cfg_p = &cfg_double[0];
		Map<VectorXd> cfg_vec(cfg_p, cfg_double.size());
		return cfg_vec;
	}
	// - BECR::LadderFact calculates the factorial contribution of the ratio from total boson number.
	double BECR::LadderFact(int n_new, int n_old)
	{
		double Factor = 1;
		if (n_new < n_old)
		{
			for (int n = n_old; n > n_new; n--)
			{
				Factor *= sqrt(n);
			}
		}
		else if (n_new > n_old)
		{
			for (int n = n_new; n > n_old; n--)
			{
				Factor /= sqrt(n);
			}
		}
		return Factor;
	}
	// - BECR::PsiUpdate will update the variational parameters of the wavefunction when given a vector of changes.
	void BECR::PsiUpdate(VectorXd dP)
	{
		return;
	}
	// - BECR::PrepPsi will load any local configuration information used in the wavefunction constituents.
	void BECR::PrepPsi(Config* Cfg)
	{
		VectorXd cfg_vec = CfgRead(Cfg);
		Occ = cfg_vec;
		OccP = Occ;
		Nb = Cfg->Nparticle();
		NbP = Nb;
		return;
	}
	// - BECR::PsiCfgUpdate will update any local configuration information after a configuration change.
	void BECR::PsiCfgUpdate() // Proposed local information changes contained within each Modifier.
	{
		Occ = OccP; Nb = NbP;
		return;
	}
	// If required later, will also have PsiCfgUpdateExt which uses an external local information source.
	// - BECR::PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
	// -- appropriate local information for one configuration loaded into the Ansatz.
	double BECR::PsiRatio(Diff diff)
	{
		double Ratio = 1;
		OccP = Occ;
		NbP = Nb;
		for (int d = 0; d < diff.num; d++)
		{
			NbP += diff.val[d];
			OccP(diff.pos[d]) += diff.val[d];
			Ratio *= pow(SPO(diff.pos[d]), diff.val[d]);
			Ratio *= LadderFact((int)OccP(diff.pos[d]), (int)Occ(diff.pos[d]));
		}
		Ratio *= LadderFact(NbP, Nb);
		return Ratio;
	}
	// If required later, will also have PsiRatioExt which loads local information to an external object.
	// - BECR::LogDeriv will return the logarithmic derivatives of the wavefunction w.r.t. its parameters.
	VectorXd BECR::LogDeriv(Config* Cfg) const
	{
		VectorXd dLogp(1);
		dLogp(0) = 0;
		return dLogp;
	}
}

#endif