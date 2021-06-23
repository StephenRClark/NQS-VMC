// ------------------------------------------------------------------------------------------------------------
// Modifier.cpp - contains the functions of the Modifier object subclasses. Modifiers apply amplitude
//				  modulations to the wavefunction that are diagonal in the configuration basis. Multiple 
//				  Modifiers can be applied to the same Ansatz.
// ------------------------------------------------------------------------------------------------------------

#ifndef MODIFIER_CPP
#define MODIFIER_CPP

#include "Modifier.hh"

namespace nqsvmc
{
	// <<<<<<------------ Gutz functions ------------>>>>>>
	// Constructors for the bosonic Gutzwiller factor:
		// - Initialise with the input Gutzwiller factor gfac.
	Gutz::Gutz(Hilbert* hlb, double gfac, bool VarFlag)
	{
		if ((h_->Type() == 'f') || (h_->Type() == 's'))
		{
			cerr << "Gutz is not currently compatible with fermionic or spin systems." << endl;
			std::abort();
		}
		N = h_->SysSize();
		Nmean = 1; // Recalculated in PrepPsi.
		VectorXd gv = VectorXd::Ones(1);
		gv(0) = gfac;
		ParamLoad(gv);
		VFlag = VarFlag;
		OptInds(0) = VarFlag;
		return;
	}
	// Internal parameter organisation functions.
	// - Gutz::ParamCheck will test the values of all parameters and ensure that they are less than the cap.
	void Gutz::ParamCheck()
	{
		if (abs(G) > ParamCap)
		{
			G = G * ParamCap / abs(G);
		}
		if (isnan(G) || isinf(G))
		{
			G = 0;
		}
		return;
	}
	// Observer functions:
	// - Gutz::Nsite returns the number of visible sites.
	int Gutz::Nsite() const
	{
		return N;
	}
	// - Gutz::Nparam returns the number of parameters.
	int Gutz::Nparam() const
	{
		return Np;
	}
	// - Gutz::VarFlag will return 0 or 1 according to whether variational modification is permitted.
	bool Gutz::VarFlag() const
	{
		return VFlag;
	}
	// - Gutz::OptIndRead returns the individual parameter flag vector.
	VectorXd Gutz::OptIndRead() const
	{
		return OptInds;
	}
	// - Gutz::ParamLoad will change the existing parameters to the provided parameters.
	void Gutz::ParamLoad(VectorXd g_new)
	{
		G = g_new(0);
		ParamCheck();
		return;
	}
	// - Gutz::ParamList returns the parameters in the Modifier
	// -- List order is <G>.
	VectorXd Gutz::ParamList() const
	{
		VectorXd Params(1);
		Params(0) = G;
		return Params;
	}

	// Variational modification management functions:
	// - Gutz::VarSwitch will alter VFlag if allowed.
	void Gutz::VarSwitch()
	{
		VFlag = !VFlag;
		return;
	}
	// - Gutz::OptIndLoad will load a Boolean vector of updated parameter flags.
	void Gutz::OptIndLoad(VectorXd newinds)
	{
		if (newinds.size() != Np)
		{
			cerr << "New optimisation index vector does not have the correct number of entries." << endl;
			std::abort();
		}
		OptInds = newinds;
	}
	// - Gutz::RndBtcSelect will randomly alter the individual parameter flags.
	// -- Default zero argument version will select 1/e of the parameters.
	void Gutz::RndBtcSelect()
	{
		return; // Only one parameter, and a minimum of one is needed.
	}
	// -- Second version with fraction will select input fraction if valid.
	void Gutz::RndBtcSelect(double pfrac)
	{
		return; // Only one parameter, and a minimum of one is needed.
	}
	// - Gutz::SetParamCap will change the maximum allowed amplitude of a single parameter.
	void Gutz::SetParamCap(double newcap)
	{
		if (newcap <= 0)
		{
			cerr << "Parameter magnitude cap must be positive definite." << endl;
			std::abort();
		}
		ParamCap = newcap;
		ParamCheck();
		return;
	}

	// Wavefunction calculation functions:
	// - Gutz::CfgRead will allow for additional modification of FullCfg methods in Config.
	VectorXd Gutz::CfgRead(Config* Cfg) const
	{
		vector<int> cfg_init = Cfg->FullCfg();
		vector<double> cfg_double(cfg_init.begin(), cfg_init.end());
		double* cfg_p = &cfg_double[0];
		Map<VectorXd> cfg_vec(cfg_p, cfg_double.size());
		return cfg_vec;
	}
	// - Gutz::PsiUpdate will update the variational parameters of the wavefunction when given a vector of changes.
	void Gutz::PsiUpdate(VectorXd dP)
	{
		G += dP(0);
		ParamCheck();
		return;
	}
	// - Gutz::PrepPsi will load any local configuration information used in the wavefunction constituents.
	void Gutz::PrepPsi(Config* Cfg)
	{
		VectorXd cfg_vec = CfgRead(Cfg);
		Nmean = cfg_vec.sum() / N;
		dN = cfg_vec.array() - Nmean;
		dNP = dN;
		return;
	}
	// - Gutz::PsiCfgUpdate will update any local configuration information after a configuration change.
	void Gutz::PsiCfgUpdate() // This version uses the stored alternate local information.
	{
		dN = dNP;
		return;
	}
	// - Gutz::PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
	// -- appropriate local information for one configuration is loaded into the Ansatz.
	double Gutz::PsiRatio(Diff diff)
	{
		dNP = dN;
		for (int d = 0; d < diff.num; d++)
		{
			dNP(diff.pos[d]) += diff.val[d];
		}
		VectorXd dN2 = (dNP.array().square() - dN.array().square());
		double Ratio = exp(-G * dN2.sum());
		if (isnan(Ratio) || isinf(Ratio))
		{
			Ratio = 0;
		}
		return Ratio;
	}
	// - Gutz::LogDeriv will return the logarithmic derivatives of the wavefunction w.r.t. its parameters.
	VectorXd Gutz::LogDeriv(Config* Cfg) const
	{
		VectorXd dLogp = VectorXd::Zero(1);
		VectorXd cfg_vec = CfgRead(Cfg);
		VectorXd dNLoc = cfg_vec.array() - Nmean;
		dNLoc = dNLoc.array().square();
		if (OptInds(0))
		{
			dLogp(0) = dNLoc.sum();
		}
		if (isnan(dLogp(0)) || isinf(dLogp(0)))
		{
			dLogp(0) = 0;
		}
		return dLogp;
	}

	// <<<<<<------------ Jast functions ------------>>>>>>
	// Constructors for the bosonic Jastrow state:
		// - Random initialisation with starting values of each parameter type.
	Jast::Jast(Hilbert* hlb, Graph* grp, double J0, vector<double> noise)
	{
		MatrixInit(hlb, grp);
		double nmag = noise[0]; // Noise vector should be populated in order: nmag, nphs.
		random_device rd;
		default_random_engine rndgen(rd());
		uniform_real_distribution<double> noisedist(-nmag, nmag);
		// Initialise random parameters, no random phase initially.
		for (int j = 0; j < Np; j++)
		{
			JsVar(j) = J0 * (1 + noisedist(rndgen));
		}
		// Starting parameter cap of 5.
		ParamCap = 5;
		ParamFill();
		// Variational flag automatically set to 1.
		VFlag = 1;
		// Automatically populate OptInds with ones.
		OptInds.resize(Np);
		OptInds.setOnes();
		return;
	}
	// - Initialisation with pre-existing parameters.
	Jast::Jast(Hilbert* hlb, Graph* grp, VectorXd ParameterList)
	{
		MatrixInit(hlb, grp);
		// Starting parameter cap of 5, or maximum given value.
		double MaxParam = max(abs(ParameterList.maxCoeff()), abs(ParameterList.minCoeff()));
		// Variational flag automatically set to 1.
		VFlag = 1;
		// Automatically populate OptInds with ones.
		OptInds.resize(Np);
		OptInds.setOnes();
		ParamCap = max(5.0, MaxParam);
		ParamLoad(ParameterList);
		return;
	}
	// Internal parameter organisation functions.
	// - Jast::MatrixInit will initialise the Jastrow parameter index matrix.
	void Jast::MatrixInit(Hilbert* hlb, Graph* grp)
	{
		h_ = hlb;
		g_ = grp;
		if ((h_->Type() == 'f') || (h_->Type() == 's'))
		{
			cerr << "This Jastrow implementation is not currently compatible with fermionic or spin systems." << endl;
			std::abort();
		}
		if (h_->SysSize() != g_->Nsite())
		{
			cerr << "Hilbert and Graph lattice size mismatch." << endl;
			std::abort();
		}
		N = h_->SysSize();
		Js = MatrixXd::Zero(N, N);
		JsI = -MatrixXd::Ones(N, N); // Use -1 to indicate unfilled spots.
		int Ntr = g_->Ntranslate();
		int pcount = 0;
		vector<int> Bond;
		int Ind1, Ind2;
		for (int n = 0; n < N; n++)
		{
			for (int m = n; m < N; m++)
			{
				if (JsI(n, m) < 0)
				{

					JsI(n, m) = pcount;
					JsI(m, n) = pcount;
					for (int b = 0; b < Ntr; b++)
					{
						Bond = g_->BondSearch(b);
						Ind1 = Bond[n];
						Ind2 = Bond[m];
						if ((Ind1 >= 0) && (Ind2 >= 0))
						{
							JsI(Ind1, Ind2) = pcount;
							JsI(Ind2, Ind1) = pcount;
						}
					}
					pcount++;
				}
			}
		}
		Np = pcount;
		JsVar.resize(Np);
		Tj.resize(N);
		TjP.resize(N);
		return;
	}
	// - Jast::ParamFill will populate the vectors and matrices given relevant parameters.
	void Jast::ParamFill()
	{
		int ind = 0;
		for (int n = 0; n < N; n++)
		{
			for (int m = 0; m < N; m++)
			{
				if (JsI(n, m) >= 0)
				{
					ind = (int)JsI(n, m);
					Js(n, m) = JsVar(ind);
				}
			}
		}
		return;
	}
	// - Jast::ParamCheck will test the values of all parameters and ensure that they are less than the cap.
	void Jast::ParamCheck()
	{
		for (int j = 0; j < Np; j++)
		{
			if (abs(JsVar(j)) > ParamCap)
			{
				JsVar(j) = JsVar(j) * ParamCap / abs(JsVar(j));
			}
			if (isnan(JsVar(j)) || isinf(JsVar(j)))
			{
				JsVar(j) = 0;
			}
		}
		return;
	}
	// Observer functions:
	// - Jast::Nsite returns the number of visible sites.
	int Jast::Nsite() const
	{
		return N;
	}
	// - Jast::Nparam returns the number of parameters.
	int Jast::Nparam() const
	{
		return Np;
	}
	// - Jast::VarFlag will return 0 or 1 according to whether variational modification is permitted.
	bool Jast::VarFlag() const
	{
		return VFlag;
	}
	// - Jast::OptIndRead returns the individual parameter flag vector.
	VectorXd Jast::OptIndRead() const
	{
		return OptInds;
	}
	// - Jast::ParamLoad will change the existing parameters to the provided parameters.
	void Jast::ParamLoad(VectorXd NewParams)
	{
		if (NewParams.size() != Np)
		{
			cerr << "Number of provided parameters (" << NewParams.size() << ") does not match Modifier parameter number (" << Np << ")." << endl;
			std::abort();
		}
		for (int p = 0; p < Np; p++)
		{
			if (abs(NewParams(p)) < 1e-30)
			{
				NewParams(p) = 0;
			}
		}
		JsVar = NewParams;
		double MaxParam = max(abs(NewParams.maxCoeff()), abs(NewParams.minCoeff()));
		ParamCap = max(ParamCap, MaxParam);
		ParamCheck();
		ParamFill();
		return;
	}
	// - Jast::ParamList returns the parameters in the Modifier.
	VectorXd Jast::ParamList() const
	{
		return JsVar;
	}

	// Variational modification management functions:
	// - Jast::VarSwitch will alter VFlag if allowed.
	void Jast::VarSwitch()
	{
		VFlag = !VFlag;
		return;
	}
	// - Jast::OptIndLoad will load a Boolean vector of updated parameter flags.
	void Jast::OptIndLoad(VectorXd newinds)
	{
		if (newinds.size() != Np)
		{
			cerr << "New optimisation index vector does not have the correct number of entries." << endl;
			std::abort();
		}
		OptInds = newinds;
		return;
	}
	// - Jast::RndBtcSelect will randomly alter the individual parameter flags.
	// -- Default zero argument version will select 1/e of the parameters.
	void Jast::RndBtcSelect()
	{
		double frac = exp(-1);
		int Num = (int)ceil(Np * frac);
		uniform_int_distribution<int> Pdist(Np);
		VectorXd newinds(Np);
		random_device rd;
		std::default_random_engine rnd(rd());
		int site;
		while (Num > 0)
		{
			site = Pdist(rnd);
			if (!newinds(site))
			{
				newinds(site) = 1;
				Num -= 1;
			}
		}
		OptIndLoad(newinds);
		return;
	}
	// -- Second version with fraction will select input fraction if valid.
	void Jast::RndBtcSelect(double pfrac)
	{
		if ((pfrac <= 0) || (pfrac >= 1))
		{
			cerr << "Invalid parameter fraction - should be between 0 and 1 non-inclusive." << endl;
			std::abort();
		}
		int Num = (int)ceil(Np * pfrac);
		uniform_int_distribution<int> Pdist(Np);
		VectorXd newinds(Np);
		random_device rd;
		std::default_random_engine rnd(rd());
		int site;
		while (Num > 0)
		{
			site = Pdist(rnd);
			if (!newinds(site))
			{
				newinds(site) = 1;
				Num -= 1;
			}
		}
		OptIndLoad(newinds);
		return;
	}
	// - Jast::SetParamCap will change the maximum allowed amplitude of a single parameter.
	void Jast::SetParamCap(double newcap)
	{
		if (newcap <= 0)
		{
			cerr << "Parameter magnitude cap must be positive definite." << endl;
			std::abort();
		}
		ParamCap = newcap;
		ParamCheck();
		ParamFill();
		return;
	}

	// Wavefunction calculation functions:
	// - Jast::CfgRead will allow for additional modification of FullCfg methods in Config.
	VectorXd Jast::CfgRead(Config* Cfg) const
	{
		vector<int> cfg_init = Cfg->FullCfg();
		vector<double> cfg_double(cfg_init.begin(), cfg_init.end());
		double* cfg_p = &cfg_double[0];
		Map<VectorXd> cfg_vec(cfg_p, cfg_double.size());
		return cfg_vec;
	}
	// - Jast::PsiUpdate will update the variational parameters of the wavefunction when given a vector of changes.
	void Jast::PsiUpdate(VectorXd dP)
	{
		dP = dP.array() * OptInds.array();
		for (int p = 0; p < Np; p++)
		{
			if (abs(dP(p)) < 1e-30)
			{
				dP(p) = 0;
			}
		}
		JsVar += dP;
		ParamCheck();
		ParamFill();
		return;
	}
	// - Jast::PrepPsi will load any local configuration information used in the wavefunction constituents.
	void Jast::PrepPsi(Config* Cfg)
	{
		VectorXd cfg_vec = CfgRead(Cfg);
		Tj = Js * cfg_vec;
		TjP = Tj;
		return;
	}
	// - Jast::PsiCfgUpdate will update any local configuration information after a configuration change.
	void Jast::PsiCfgUpdate() // This version uses the stored alternate local information.
	{
		Tj = TjP;
		return;
	}
	// - Jast::PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
	// -- appropriate local information for one configuration is loaded into the Ansatz.
	double Jast::PsiRatio(Diff diff)
	{
		double Ratio = 1;
		TjP = Tj;
		VectorXd DeltaVec = VectorXd::Zero(N);
		for (int d = 0; d < diff.num; d++)
		{
			DeltaVec(diff.pos[d]) = diff.val[d];
			Ratio *= exp(-Tj(diff.pos[d]) * diff.val[d]);
			for (int n = 0; n < N; n++)
			{
				TjP(n) += diff.val[d] * Js(diff.pos[d], n);
			}
		}
		MatrixXd DeltaMat = DeltaVec * DeltaVec.transpose();
		DeltaMat = DeltaMat.array() * Js.array();
		Ratio *= exp(-DeltaMat.sum() / 2);
		if (isnan(Ratio) || isinf(Ratio))
		{
			Ratio = 0;
		}
		return Ratio;
	}
	// - Jast::LogDeriv will return the logarithmic derivatives of the wavefunction w.r.t. its parameters.
	VectorXd Jast::LogDeriv(Config* Cfg) const
	{
		VectorXd dLogp = VectorXd::Zero(Np);
		VectorXd cfg_vec = CfgRead(Cfg);
		MatrixXd cfg_mat = cfg_vec * cfg_vec.transpose();
		for (int n = 0; n < N; n++)
		{
			for (int m = 0; m < N; m++)
			{
				if (JsI(n, m) >= 0)
				{
					dLogp((int)JsI(n, m)) -= 0.5 * cfg_mat(n, m);
				}
			}
		}
		return dLogp;
	}

	// <<<<<<------------ MBC functions ------------>>>>>>
	// Constructors for the bosonic many body correlator:
	// - Random initialisation with starting values of each parameter type.
	NNMB::NNMB(Hilbert* hlb, Graph* grp, double G0, vector<double> noise)
	{
		ListInit(hlb, grp);
		double nmag = noise[0]; // Noise vector should be populated in order: nmag, nphs.
		random_device rd;
		default_random_engine rndgen(rd());
		uniform_real_distribution<double> noisedist(-nmag, nmag);
		// Initialise random parameters, no random phase initially.
		gmb = G0 * (1 + noisedist(rndgen));
		// Starting parameter cap of 5.
		ParamCap = 5;
		ParamCheck();
		// Variational flag automatically set to 1.
		VFlag = 1;
		// Automatically populate OptInds with ones.
		OptInds.resize(Np);
		OptInds.setOnes();
		return;
	}
	// - Initialisation with pre-existing parameters.
	NNMB::NNMB(Hilbert* hlb, Graph* grp, VectorXd ParameterList)
	{
		ListInit(hlb, grp);
		// Starting parameter cap of 5, or maximum given value.
		double MaxParam = max(abs(ParameterList.maxCoeff()), abs(ParameterList.minCoeff()));
		// Variational flag automatically set to 1.
		VFlag = 1;
		// Automatically populate OptInds with ones.
		OptInds.resize(Np);
		OptInds.setOnes();
		ParamCap = max(5.0, MaxParam);
		ParamLoad(ParameterList);
		return;
	}
	// Internal parameter organisation functions.
	// - NNMB::ListInit will initialise the nearest neighbour index lists.
	void NNMB::ListInit(Hilbert* hlb, Graph* grp)
	{
		h_ = hlb;
		g_ = grp;
		Np = 1;
		if ((h_->Type() == 'f') || (h_->Type() == 's'))
		{
			cerr << "NNMB is not currently compatible with fermionic or spin systems." << endl;
			std::abort();
		}
		if (h_->SysSize() != g_->Nsite())
		{
			cerr << "Hilbert and Graph lattice size mismatch." << endl;
			std::abort();
		}
		N = h_->SysSize();
		z = g_->Coord();
		int nvec = g_->Nvecs();
		vector<vector<int>> Bonds = g_->BondRead();
		NNList.resize(N);
		for (int n = 0; n < N; n++)
		{
			NNList[n].resize(z);
			for (int b = 0; b < nvec; b++)
			{
				NNList[n][b] = Bonds[b][n];
			}
		}
		dN.resize(N);
		dNP.resize(N);
		Xi.resize(N);
		XiP.resize(N);
		return;
	}
	// - NNMB::ParamCheck will test the values of all parameters and ensure that they are less than the cap.
	void NNMB::ParamCheck()
	{
		if (abs(gmb) > ParamCap)
		{
			gmb = gmb * ParamCap / abs(gmb);
		}
		if (isnan(gmb) || isinf(gmb))
		{
			gmb = 0;
		}
		return;
	}
	// Observer functions:
	// - NNMB::Nsite returns the number of visible sites.
	int NNMB::Nsite() const
	{
		return N;
	}
	// - NNMB::Nparam returns the number of parameters.
	int NNMB::Nparam() const
	{
		return 1;
	}
	// - NNMB::VarFlag will return 0 or 1 according to whether variational modification is permitted.
	bool NNMB::VarFlag() const
	{
		return VFlag;
	}
	// - NNMB::OptIndRead returns the individual parameter flag vector.
	VectorXd NNMB::OptIndRead() const
	{
		return OptInds;
	}
	// - NNMB::ParamLoad will change the existing parameters to the provided parameters.
	void NNMB::ParamLoad(VectorXd g_new)
	{
		gmb = g_new(0);
		ParamCheck();
		return;
	}
	// - NNMB::ParamList returns the parameters in the Modifier.
	VectorXd NNMB::ParamList() const
	{
		VectorXd gv = VectorXd::Zero(1);
		gv(0) = gmb;
		return gv;
	}

	// Variational modification management functions:
	// - NNMB::VarSwitch will alter VFlag if allowed.
	void NNMB::VarSwitch()
	{
		VFlag = !VFlag;
		return;
	}
	// - NNMB::OptIndLoad will load a Boolean vector of updated parameter flags.
	void NNMB::OptIndLoad(VectorXd newinds)
	{
		if (newinds.size() != Np)
		{
			cerr << "New optimisation index vector does not have the correct number of entries." << endl;
			std::abort();
		}
		OptInds = newinds;
	}
	// - NNMB::RndBtcSelect will randomly alter the individual parameter flags.
	// -- Default zero argument version will select 1/e of the parameters.
	void NNMB::RndBtcSelect()
	{
		return; // Only one parameter, and a minimum of one is needed.
	}
	// -- Second version with fraction will select input fraction if valid.
	void NNMB::RndBtcSelect(double pfrac)
	{
		return; // Only one parameter, and a minimum of one is needed.
	}
	// - NNMB::SetParamCap will change the maximum allowed amplitude of a single parameter.
	void NNMB::SetParamCap(double newcap)
	{
		if (newcap <= 0)
		{
			cerr << "Parameter magnitude cap must be positive definite." << endl;
			std::abort();
		}
		ParamCap = newcap;
		ParamCheck();
		return;
	}

	// Wavefunction calculation functions:
	// - NNMB::CfgRead will allow for additional modification of FullCfg methods in Config.
	VectorXd NNMB::CfgRead(Config* Cfg) const
	{
		vector<int> cfg_init = Cfg->FullCfg();
		vector<double> cfg_double(cfg_init.begin(), cfg_init.end());
		double* cfg_p = &cfg_double[0];
		Map<VectorXd> cfg_vec(cfg_p, cfg_double.size());
		return cfg_vec;
	}
	// - NNMB::PsiUpdate will update the variational parameters of the wavefunction when given a vector of changes.
	void NNMB::PsiUpdate(VectorXd dP)
	{
		gmb += dP(0);
		ParamCheck();
		return;
	}
	// - NNMB::PrepPsi will load any local configuration information used in the wavefunction constituents.
	void NNMB::PrepPsi(Config* Cfg)
	{
		VectorXd cfg_vec = CfgRead(Cfg);
		Nmean = round((cfg_vec.sum() / N));
		dN = cfg_vec.array() - Nmean;
		for (int n = 0; n < N; n++)
		{
			if ((dN(n) == 1) || (dN(n) == -1))
			{
				Xi(n) = 1;
				for (int t = 0; t < z; t++)
				{
					int nn = NNList[n][t];
					if (nn >= 0)
					{
						if (dN(nn) == -dN(nn))
						{
							Xi(n) = 0;
						}
					}
				}
			}
		}
		dNP = dN;
		XiP = Xi;
		return;
	}
	// - NNMB::PsiCfgUpdate will update any local configuration information after a configuration change.
	void NNMB::PsiCfgUpdate() // This version uses the stored alternate local information.
	{
		Xi = XiP; dN = dNP;
		return;
	}
	// - NNMB::PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
	// -- appropriate local information for one configuration is loaded into the Ansatz.
	double NNMB::PsiRatio(Diff diff)
	{
		dNP = dN; XiP = Xi;
		VectorXd EList = VectorXd::Zero(N); // List of sites requiring operator recalculation.
		for (int d = 0; d < diff.num; d++)
		{
			dNP(diff.pos[d]) += diff.val[d];
			EList(diff.pos[d]) = 1;
			for (int t = 0; t < z; t++)
			{
				int nn = NNList[diff.pos[d]][t];
				if (nn >= 0)
				{
					EList(nn) = 1;
				}
			}
		}
		for (int n = 0; n < N; n++)
		{
			if (EList(n) == 1)
			{
				if ((dNP(n) == 1) || (dNP(n) == -1))
				{
					XiP(n) = 1;
					for (int nn = 0; nn < z; nn++)
					{
						int m = NNList[n][nn];
						if (m >= 0)
						{
							if (dNP(m) == -dNP(n))
							{
								XiP(n) = 0;
							}
						}
					}
				}

			}
		}
		VectorXd dXi = XiP - Xi;
		double Ratio = exp(gmb * dXi.sum());
		if (isnan(Ratio) || isinf(Ratio))
		{
			Ratio = 0;
		}
		return Ratio;
	}
	// - NNMB::LogDeriv will return the logarithmic derivatives of the wavefunction w.r.t. its parameters.
	VectorXd NNMB::LogDeriv(Config* Cfg) const
	{
		VectorXd dLogp = VectorXd::Zero(Np);
		if (OptInds(0))
		{
			dLogp(0) = Xi.sum();
		}
		if (isnan(dLogp(0)) || isinf(dLogp(0)))
		{
			dLogp(0) = 0;
		}
		return dLogp;
	}
}

#endif