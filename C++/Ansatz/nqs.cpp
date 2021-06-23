// ------------------------------------------------------------------------------------------------------------
// nqs.cpp - contains the functions of the NQS Modifier subclasses.
// ------------------------------------------------------------------------------------------------------------

#ifndef NQS_CPP
#define NQS_CPP

#include "nqs.hh"

namespace nqsvmc
{
	// <<<<<<------------ NQSSH functions ------------>>>>>>
	// Constructors for the spin-hidden NQS:
	// - Random initialisation with starting values of each parameter type.
	NQSSH::NQSSH(Hilbert* hlb, Graph* grp, int HUDen, vector<double> StartParams, vector<double> noise)
	{
		NetworkInit(hlb, grp, HUDen);
		double nmag = noise[0]; // Noise vector should be populated in order: nmag, nphs.
		random_device rd;
		default_random_engine rndgen(rd());
		uniform_real_distribution<double> noisedist(-nmag, nmag);
		// Initialise random parameters, no random phase initially.
		a = StartParams[0] * (1 + noisedist(rndgen));
		A = StartParams[1] * (1 + noisedist(rndgen));
		for (int a = 0; a < Alpha; a++)
		{
			b(a) = StartParams[2] * (1 + noisedist(rndgen));
			B(a) = StartParams[3] * (1 + noisedist(rndgen));
			for (int n = 0; n < Nv; n++)
			{
				Wv(a, n) = StartParams[4] * (1 + noisedist(rndgen));
			}
		}
		// Automatically populate OptInds with ones.
		OptInds.resize(Np);
		OptInds.setOnes();
		// Starting parameter cap of 5.
		ParamCap = 5;
		ParamFill();
		// Variational flag automatically set to 1.
		VFlag = 1;
		return;
	}
	// - Initialisation with pre-existing parameters.
	NQSSH::NQSSH(Hilbert* hlb, Graph* grp, int HUDen, VectorXd ParameterList)
	{
		NetworkInit(hlb, grp, HUDen);
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
	// - NQSSH::NetworkInit will initialise the visible and hidden units of the RBM.
	void NQSSH::NetworkInit(Hilbert* hlb, Graph* grp, int HUDen)
	{
		h_ = hlb;
		g_ = grp;
		if (h_->Type() == 'f')
		{
			cerr << "NQSSH is not compatible with fermionic systems." << endl;
			std::abort();
		}
		if (h_->SysSize() != g_->Nsite())
		{
			cerr << "Hilbert and Graph lattice size mismatch." << endl;
			std::abort();
		}
		Alpha = HUDen;
		HDim = h_->SiteDim();
		Nv = h_->SysSize();
		Nh = Alpha * g_->Ntranslate();
		Np = 2 + 2 * Alpha + Alpha * Nv;
		av.resize(Nv);
		Av.resize(Nv);
		b.resize(Alpha);
		B.resize(Alpha);
		bv.resize(Nh);
		Bv.resize(Nh);
		Wv.resize(Alpha, Nv);
		Wm.resize(Nh, Nv);
		Theta.resize(Nh);
		Nsq.resize(Nv);
		return;
	}
	// - NQSSH::ParamFill will populate the vectors and matrices given relevant parameters.
	void NQSSH::ParamFill()
	{
		int Ntr = g_->Ntranslate();
		for (int n = 0; n < Nv; n++)
		{
			av(n) = a;
			Av(n) = A;
		}
		for (int h = 0; h < Alpha; h++)
		{
			for (int t = 0; t < Ntr; t++)
			{
				int HInd = t + h * Ntr;
				bv(HInd) = b(h);
				Bv(HInd) = B(h);
				vector<int> Bond = g_->BondSearch(t);
				for (int n = 0; n < Nv; n++)
				{
					if (Bond[n] >= 0)
					{
						Wm(HInd, Bond[n]) = Wv(h, n);
					}
				}
			}
		}
		return;
	}
	// - NQSSH::ParamCheck will test the values of all parameters and ensure that they are less than the cap.
	void NQSSH::ParamCheck()
	{
		if (abs(a) > ParamCap)
		{
			a = a * ParamCap / abs(a);
		}
		if (abs(A) > ParamCap)
		{
			A = A * ParamCap / abs(A);
		}
		if (isnan(a) || isinf(a))
		{
			a = 0;
		}
		if (isnan(A) || isinf(A))
		{
			A = 0;
		}
		for (int al = 0; al < Alpha; al++)
		{
			if (abs(b(al)) > ParamCap)
			{
				b(al) = b(al) * ParamCap / abs(b(al));
			}
			if (abs(B(al)) > ParamCap)
			{
				B(al) = B(al) * ParamCap / abs(B(al));
			}
			if (isnan(b(al)) || isinf(b(al)))
			{
				b(al) = 0;
			}
			if (isnan(B(al)) || isinf(B(al)))
			{
				B(al) = 0;
			}
			for (int n = 0; n < Nv; n++)
			{
				if (abs(Wv(al, n)) > ParamCap)
				{
					Wv(al, n) = Wv(al, n) * ParamCap / abs(Wv(al, n));
				}
				if (isnan(Wv(al, n)) || isinf(Wv(al, n)))
				{
					Wv(al, n) = 0;
				}
			}
		}
		return;
	}
	// Observer functions:
	// - NQSSH::Nvisible returns the number of visible sites.
	int NQSSH::Nvisible() const
	{
		return Nv;
	}
	// - NQSSH::Nhidden returns the number of hidden sites.
	int NQSSH::Nhidden() const
	{
		return Nh;
	}
	// - NQSSH::Nparam returns the number of parameters.
	int NQSSH::Nparam() const
	{
		return Np;
	}
	// - NQSSH::VarFlag will return 0 or 1 according to whether variational modification is permitted.
	bool NQSSH::VarFlag() const
	{
		return VFlag;
	}
	// - NQSSH::OptIndRead returns the individual parameter flag vector.
	VectorXd NQSSH::OptIndRead() const
	{
		return OptInds;
	}
	// - NQSSH::ParamLoad will change the existing parameters to the provided parameters.
	void NQSSH::ParamLoad(VectorXd NewParams)
	{
		if (NewParams.size() != Np)
		{
			cerr << "Number of provided parameters does not match Modifier parameter number." << endl;
			std::abort();
		}
		for (int p = 0; p < Np; p++)
		{
			if (abs(NewParams(p)) < 1e-30)
			{
				NewParams(p) = 0;
			}
		}
		NewParams = NewParams.array() * OptInds.array();
		double da = NewParams(0);
		double dA = NewParams(1);
		VectorXd db = NewParams.segment(2, Alpha);
		VectorXd dB = NewParams.segment(2 + Alpha, Alpha);
		VectorXd dW = NewParams.tail(Alpha * Nv);
		MatrixXd dWv(Map<MatrixXd>(dW.data(), Nv, Alpha));
		dWv.transposeInPlace();
		a = da;
		A = dA;
		b = db;
		B = dB;
		Wv = dWv;
		double MaxParam = max(abs(NewParams.maxCoeff()), abs(NewParams.minCoeff()));
		ParamCap = max(ParamCap, MaxParam);
		ParamCheck();
		ParamFill();
		return;
	}
	// - NQSSH::ParamList returns the parameters in the Modifier
	// -- List order is <a, A, b, B, Wv>.
	VectorXd NQSSH::ParamList() const
	{
		VectorXd Params;
		Params.resize(Np);
		Params(0) = a;
		Params(1) = A;
		Params.segment(2, Alpha) = b;
		Params.segment(2 + Alpha, Alpha) = B;
		MatrixXd Wt = Wv.transpose();
		VectorXd Wvec(Map<VectorXd>(Wt.data(), Wv.cols() * Wv.rows()));
		Params.tail(Alpha * Nv) = Wvec;
		return Params;
	}

	// Variational modification management functions:
	// - NQSSH::VarSwitch will alter VFlag if allowed.
	void NQSSH::VarSwitch()
	{
		VFlag = !VFlag;
		return;
	}
	// - NQSSH::OptIndLoad will load a Boolean vector of updated parameter flags.
	void NQSSH::OptIndLoad(VectorXd newinds)
	{
		if (newinds.size() != Np)
		{
			cerr << "New optimisation index vector does not have the correct number of entries." << endl;
			std::abort();
		}
		OptInds = newinds;
		return;
	}
	// - NQSSH::RndBtcSelect will randomly alter the individual parameter flags.
	// -- Default zero argument version will select 1/e of the parameters.
	void NQSSH::RndBtcSelect()
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
	void NQSSH::RndBtcSelect(double pfrac)
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
	// - NQSSH::SetParamCap will change the maximum allowed amplitude of a single parameter.
	void NQSSH::SetParamCap(double newcap)
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
	VectorXd NQSSH::SHTrace(VectorXd theta_in) const // Calculate trace over hidden units for given Theta.
	{
		VectorXd h = VectorXd::LinSpaced(HDim, (1 - HDim) / 2, (HDim - 1) / 2);
		VectorXd hsq = h.array().square();
		MatrixXd F = theta_in * h.transpose();
		MatrixXd Bh = Bv * hsq.transpose();
		F = F.array() + Bh.array();
		F = F.array().exp();
		VectorXd Tr = VectorXd::Zero(F.rows()); //  F.rowwise().sum() / HDim;
		for (int c = 0; c < F.cols(); c++)
		{
			Tr += (F.col(c) / (double)HDim);
		}
		return Tr;
	}
	VectorXd NQSSH::dT_SHTrace(VectorXd theta_in) const// Calculate derivative of trace w.r.t Theta.
	{
		VectorXd h = VectorXd::LinSpaced(HDim, (1 - HDim) / 2, (HDim - 1) / 2);
		VectorXd hsq = h.array().square();
		MatrixXd hm = VectorXd::Ones(theta_in.size()) * h.transpose();
		MatrixXd F = theta_in * h.transpose();
		MatrixXd Bh = Bv * hsq.transpose();
		F = F.array() + Bh.array();
		F = F.array().exp();
		MatrixXd dF = F.array() * hm.array();
		VectorXd dFv = VectorXd::Zero(F.rows());
		VectorXd Fv = VectorXd::Zero(F.rows());
		for (int c = 0; c < F.cols(); c++)
		{
			dFv += dF.col(c);
			Fv += F.col(c);
		}
		VectorXd dTr = dFv.array() / Fv.array();
		return dTr;
	}
	VectorXd NQSSH::dB_SHTrace(VectorXd theta_in) const // Calculate derivative of trace w.r.t B.
	{
		VectorXd h = VectorXd::LinSpaced(HDim, (1 - HDim) / 2, (HDim - 1) / 2);
		VectorXd hsq = h.array().square();
		MatrixXd h2m = VectorXd::Ones(theta_in.size()) * hsq.transpose();
		MatrixXd F = theta_in * h.transpose();
		MatrixXd Bh = Bv * hsq.transpose();
		F = F.array() + Bh.array();
		F = F.array().exp();
		MatrixXd dF = F.array() * h2m.array();
		VectorXd dFv = VectorXd::Zero(F.rows());
		VectorXd Fv = VectorXd::Zero(F.rows());
		for (int c = 0; c < F.cols(); c++)
		{
			dFv += dF.col(c);
			Fv += F.col(c);
		}
		VectorXd dTr = dFv.array() / Fv.array();
		return dTr;
	}
	// - NQSSH::CfgRead will allow for additional modification of FullCfg methods in Config.
	VectorXd NQSSH::CfgRead(Config* Cfg) const
	{
		vector<int> cfg_init = Cfg->FullCfg();
		vector<double> cfg_double(cfg_init.begin(), cfg_init.end());
		double* cfg_p = &cfg_double[0];
		Map<VectorXd> cfg_vec(cfg_p, cfg_double.size());
		return cfg_vec;
	}
	// - NQSSH::PsiUpdate will update the variational parameters of the wavefunction when given a vector of changes.
	void NQSSH::PsiUpdate(VectorXd dP)
	{
		dP = dP.array() * OptInds.array();
		for (int p = 0; p < Np; p++)
		{
			if (abs(dP(p)) < 1e-30)
			{
				dP(p) = 0;
			}
		}
		dP = dP.array() * OptInds.array();
		double da = dP(0);
		double dA = dP(1);
		VectorXd db = dP.segment(2, Alpha);
		VectorXd dB = dP.segment(2 + Alpha, Alpha);
		VectorXd dW = dP.tail(Alpha * Nv);
		MatrixXd dWv(Map<MatrixXd>(dW.data(), Nv, Alpha));
		dWv.transposeInPlace();
		a += da;
		A += dA;
		b += db;
		B += dB;
		Wv += dWv;
		ParamCheck();
		ParamFill();
		return;
	}
	// - NQSSH::PrepPsi will load any local configuration information used in the wavefunction constituents.
	void NQSSH::PrepPsi(Config* Cfg)
	{
		VectorXd cfg_vec = CfgRead(Cfg);
		Theta = bv + (Wm * cfg_vec);
		Nsq = cfg_vec.array().square();
		ThetaP = Theta;
		NsqP = Nsq;
		return;
	}
	// - NQSSH::PsiCfgUpdate will update any local configuration information after a configuration change.
	void NQSSH::PsiCfgUpdate() // This version uses the stored alternate local information.
	{
		Theta = ThetaP;
		Nsq = NsqP;
		return;
	}
	// - NQSSH::PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
	// -- appropriate local information for one configuration is loaded into the Ansatz.
	double NQSSH::PsiRatio(Diff diff)
	{
		double Ratio = 1;
		VectorXd ThetaShift = VectorXd::Zero(Nh); // May need to redo this initialisation.
		VectorXd NsqShift = VectorXd::Zero(Nv);
		for (int d = 0; d < diff.num; d++)
		{
			Ratio *= exp(av(diff.pos[d]) * diff.val[d]);
			ThetaShift += (Wm.col(diff.pos[d]) * diff.val[d]);
			NsqShift(diff.pos[d]) = ((2 * sqrt(Nsq(diff.pos[d])) + diff.val[d]) * diff.val[d]);
			Ratio *= exp(Av(diff.pos[d]) * NsqShift(diff.pos[d]));
		}
		NsqP = Nsq + NsqShift;
		ThetaP = Theta + ThetaShift;
		VectorXd Trace = SHTrace(Theta);
		VectorXd TraceP = SHTrace(ThetaP);
		VectorXd TraceRatio = TraceP.array() / Trace.array();
		Ratio *= TraceRatio.prod();
		if (isnan(Ratio) || isinf(Ratio))
		{
			Ratio = 0;
		}
		return Ratio;
	}
	// - NQSSH::LogDeriv will return the logarithmic derivatives of the wavefunction w.r.t. its parameters.
	VectorXd NQSSH::LogDeriv(Config* Cfg) const
	{
		VectorXd dLogp = VectorXd::Zero(Np);
		VectorXd cfg_vec = CfgRead(Cfg);
		int Ntr = g_->Ntranslate();
		if (OptInds(0))
		{
			dLogp(0) = cfg_vec.sum();
		}
		if (OptInds(1))
		{
			dLogp(1) = Nsq.sum();
		}
		VectorXd dTheta = dT_SHTrace(Theta);
		VectorXd dB = dB_SHTrace(Theta);
		for (int al = 0; al < Alpha; al++)
		{
			if (OptInds(al + 2))
			{
				VectorXd dThetaA = dTheta.segment(al * Ntr, Ntr);
				dLogp(al + 2) = dThetaA.sum();
			}
			if (OptInds(al + Alpha + 2))
			{
				VectorXd dBA = dB.segment(al * Ntr, Ntr);
				dLogp(al + 2 + Alpha) = dBA.sum();
			}
			for (int b = 0; b < Ntr; b++)
			{
				vector<int> Bond = g_->BondSearch(b);
				int TInd = b + al * Ntr;
				for (int n = 0; n < Nv; n++)
				{
					int PInd = 2 + (2 * Alpha) + n + (al * Nv);
					int VInd = Bond[n];
					if (OptInds(PInd))
					{
						dLogp(PInd) += (cfg_vec(VInd) * dTheta(TInd));
					}
				}
			}
		}
		for (int np = 0; np < Np; np++)
		{
			if (isnan(dLogp(np)) || isinf(dLogp(np)))
			{
				dLogp(np) = 0;
			}
		}
		return dLogp;
	}

	// <<<<<<------------ NQSNH functions ------------>>>>>>
	// Constructors for the number-hidden NQS:
		// - Random initialisation with starting values of each parameter type.
	NQSNH::NQSNH(Hilbert* hlb, Graph* grp, int HUDen, vector<double> StartParams, vector<double> noise)
	{
		NetworkInit(hlb, grp, HUDen);
		double nmag = noise[0]; // Noise vector should be populated in order: seed, nmag, nphs.
		random_device rd;
		default_random_engine rndgen(rd());
		uniform_real_distribution<double> noisedist(-nmag, nmag);
		// Initialise random parameters, no random phase initially.
		a = StartParams[0] * (1 + noisedist(rndgen));
		A = StartParams[1] * (1 + noisedist(rndgen));
		for (int a = 0; a < Alpha; a++)
		{
			b(a) = StartParams[2] * (1 + noisedist(rndgen));
			B(a) = StartParams[3] * (1 + noisedist(rndgen));
			for (int n = 0; n < Nv; n++)
			{
				Wv(a, n) = StartParams[4] * (1 + noisedist(rndgen));
			}
		}
		// Automatically populate OptInds with ones.
		OptInds.resize(Np);
		OptInds.setOnes();
		// Starting parameter cap of 5.
		ParamCap = 5;
		ParamFill();
		// Variational flag automatically set to 1.
		VFlag = 1;
		return;
	}
	// - Initialisation with pre-existing parameters.
	NQSNH::NQSNH(Hilbert* hlb, Graph* grp, int HUDen, VectorXd ParameterList)
	{
		NetworkInit(hlb, grp, HUDen);
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
	// - NQSNH::NetworkInit will initialise the visible and hidden units of the RBM.
	void NQSNH::NetworkInit(Hilbert* hlb, Graph* grp, int HUDen)
	{
		h_ = hlb;
		g_ = grp;
		if (h_->Type() == 'f')
		{
			cerr << "NQSSH is not compatible with fermionic systems." << endl;
			std::abort();
		}
		if (h_->SysSize() != g_->Nsite())
		{
			cerr << "Hilbert and Graph lattice size mismatch." << endl;
			std::abort();
		}
		Alpha = HUDen;
		HDim = h_->SiteDim();
		Nv = h_->SysSize();
		Nh = Alpha * g_->Ntranslate();
		Np = 2 + 2 * Alpha + Alpha * Nv;
		av.resize(Nv);
		Av.resize(Nv);
		b.resize(Alpha);
		B.resize(Alpha);
		bv.resize(Nh);
		Bv.resize(Nh);
		Wv.resize(Alpha, Nv);
		Wm.resize(Nh, Nv);
		Theta.resize(Nh);
		Nsq.resize(Nv);
		return;
	}
	// - NQSNH::ParamFill will populate the vectors and matrices given relevant parameters.
	void NQSNH::ParamFill()
	{
		int Ntr = g_->Ntranslate();
		for (int n = 0; n < Nv; n++)
		{
			av(n) = a;
			Av(n) = A;
		}
		for (int h = 0; h < Alpha; h++)
		{
			for (int t = 0; t < Ntr; t++)
			{
				int HInd = t + h * Ntr;
				bv(HInd) = b(h);
				Bv(HInd) = B(h);
				vector<int> Bond = g_->BondSearch(t);
				for (int n = 0; n < Nv; n++)
				{
					if (Bond[n] >= 0)
					{
						Wm(HInd, Bond[n]) = Wv(h, n);
					}
				}
			}
		}
		return;
	}
	// - NQSNH::ParamCheck will test the values of all parameters and ensure that they are less than the cap.
	void NQSNH::ParamCheck()
	{
		if (abs(a) > ParamCap)
		{
			a = a * ParamCap / abs(a);
		}
		if (abs(A) > ParamCap)
		{
			A = A * ParamCap / abs(A);
		}
		if (isnan(a) || isinf(a))
		{
			a = 0;
		}
		if (isnan(A) || isinf(A))
		{
			A = 0;
		}
		for (int al = 0; al < Alpha; al++)
		{
			if (abs(b(al)) > ParamCap)
			{
				b(al) = b(al) * ParamCap / abs(b(al));
			}
			if (abs(B(al)) > ParamCap)
			{
				B(al) = B(al) * ParamCap / abs(B(al));
			}
			if (isnan(b(al)) || isinf(b(al)))
			{
				b(al) = 0;
			}
			if (isnan(B(al)) || isinf(B(al)))
			{
				B(al) = 0;
			}
			for (int n = 0; n < Nv; n++)
			{
				if (abs(Wv(al, n)) > ParamCap)
				{
					Wv(al, n) = Wv(al, n) * ParamCap / abs(Wv(al, n));
				}
				if (isnan(Wv(al, n)) || isinf(Wv(al, n)))
				{
					Wv(al, n) = 0;
				}
			}
		}
		return;
	}
	// Observer functions:
	// - NQSNH::Nvisible returns the number of visible sites.
	int NQSNH::Nvisible() const
	{
		return Nv;
	}
	// - NQSNH::Nhidden returns the number of hidden sites.
	int NQSNH::Nhidden() const
	{
		return Nh;
	}
	// - NQSNH::Nparam returns the number of parameters.
	int NQSNH::Nparam() const
	{
		return Np;
	}
	// - NQSNH::VarFlag will return 0 or 1 according to whether variational modification is permitted.
	bool NQSNH::VarFlag() const
	{
		return VFlag;
	}
	// - NQSNH::OptIndRead returns the individual parameter flag vector.
	VectorXd NQSNH::OptIndRead() const
	{
		return OptInds;
	}
	// - NQSNH::ParamLoad will change the existing parameters to the provided parameters.
	void NQSNH::ParamLoad(VectorXd NewParams)
	{
		if (NewParams.size() != Np)
		{
			cerr << "Number of provided parameters does not match Modifier parameter number." << endl;
			std::abort();
		}
		for (int p = 0; p < Np; p++)
		{
			if (abs(NewParams(p)) < 1e-30)
			{
				NewParams(p) = 0;
			}
		}
		NewParams = NewParams.array() * OptInds.array();
		double da = NewParams(0);
		double dA = NewParams(1);
		VectorXd db = NewParams.segment(2, Alpha);
		VectorXd dB = NewParams.segment(2 + Alpha, Alpha);
		VectorXd dW = NewParams.tail(Alpha * Nv);
		MatrixXd dWv(Map<MatrixXd>(dW.data(), Nv, Alpha));
		dWv.transposeInPlace();
		a = da;
		A = dA;
		b = db;
		B = dB;
		Wv = dWv;
		double MaxParam = max(abs(NewParams.maxCoeff()), abs(NewParams.minCoeff()));
		ParamCap = max(ParamCap, MaxParam);
		ParamCheck();
		ParamFill();
		return;
	}

	// - ParamList returns the parameters in the Modifier
	// -- List order is <a, A, b, B, Wv>.
	VectorXd NQSNH::ParamList() const
	{
		VectorXd Params;
		Params.resize(Np);
		Params(0) = a;
		Params(1) = A;
		Params.segment(2, Alpha) = b;
		Params.segment(2 + Alpha, Alpha) = B;
		MatrixXd Wt = Wv.transpose();
		VectorXd Wvec(Map<VectorXd>(Wt.data(), Wv.cols() * Wv.rows()));
		Params.tail(Alpha * Nv) = Wvec;
		return Params;
	}

	// Variational modification management functions:
	// - NQSNH::VarSwitch will alter VFlag if allowed.
	void NQSNH::VarSwitch()
	{
		VFlag = !VFlag;
		return;
	}
	// - NQSNH::OptIndLoad will load a Boolean vector of updated parameter flags.
	void NQSNH::OptIndLoad(VectorXd newinds)
	{
		if (newinds.size() != Np)
		{
			cerr << "New optimisation index vector does not have the correct number of entries." << endl;
			std::abort();
		}
		OptInds = newinds;
		return;
	}
	// - NQSNH::RndBtcSelect will randomly alter the individual parameter flags.
	// -- Default zero argument version will select 1/e of the parameters.
	void NQSNH::RndBtcSelect()
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
	void NQSNH::RndBtcSelect(double pfrac)
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
	// - NQSNH::SetParamCap will change the maximum allowed amplitude of a single parameter.
	void NQSNH::SetParamCap(double newcap)
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
	VectorXd NQSNH::NHTrace(VectorXd theta_in) const // Calculate trace over hidden units for given Theta.
	{
		VectorXd h = VectorXd::LinSpaced(HDim, 0, (HDim - 1));
		VectorXd hsq = h.array().square();
		MatrixXd F = theta_in * h.transpose();
		MatrixXd Bh = Bv * hsq.transpose();
		F = F.array() + Bh.array();
		F = F.array().exp();
		VectorXd Tr = VectorXd::Zero(F.rows()); //  F.rowwise().sum() / HDim;
		for (int c = 0; c < F.cols(); c++)
		{
			Tr += (F.col(c) / (double)HDim);
		}
		return Tr;
	}
	VectorXd NQSNH::dT_NHTrace(VectorXd theta_in) const// Calculate derivative of trace w.r.t Theta.
	{
		VectorXd h = VectorXd::LinSpaced(HDim, 0, (HDim - 1));
		VectorXd hsq = h.array().square();
		MatrixXd hm = VectorXd::Ones(theta_in.size()) * h.transpose();
		MatrixXd F = theta_in * h.transpose();
		MatrixXd Bh = Bv * hsq.transpose();
		F = F.array() + Bh.array();
		F = F.array().exp();
		MatrixXd dF = F.array() * hm.array();
		VectorXd dFv = VectorXd::Zero(F.rows());
		VectorXd Fv = VectorXd::Zero(F.rows());
		for (int c = 0; c < F.cols(); c++)
		{
			dFv += dF.col(c);
			Fv += F.col(c);
		}
		VectorXd dTr = dFv.array() / Fv.array();
		return dTr;
	}
	VectorXd NQSNH::dB_NHTrace(VectorXd theta_in) const // Calculate derivative of trace w.r.t B.
	{
		VectorXd h = VectorXd::LinSpaced(HDim, 0, (HDim - 1));
		VectorXd hsq = h.array().square();
		MatrixXd h2m = VectorXd::Ones(theta_in.size()) * hsq.transpose();
		MatrixXd F = theta_in * h.transpose();
		MatrixXd Bh = Bv * hsq.transpose();
		F = F.array() + Bh.array();
		F = F.array().exp();
		MatrixXd dF = F.array() * h2m.array();
		VectorXd dFv = VectorXd::Zero(F.rows());
		VectorXd Fv = VectorXd::Zero(F.rows());
		for (int c = 0; c < F.cols(); c++)
		{
			dFv += dF.col(c);
			Fv += F.col(c);
		}
		VectorXd dTr = dFv.array() / Fv.array();
		return dTr;
	}
	// - NQSNH::CfgRead will allow for additional modification of FullCfg methods in Config.
	VectorXd NQSNH::CfgRead(Config* Cfg) const
	{
		vector<int> cfg_init = Cfg->FullCfg();
		vector<double> cfg_double(cfg_init.begin(), cfg_init.end());
		double* cfg_p = &cfg_double[0];
		Map<VectorXd> cfg_vec(cfg_p, cfg_double.size());
		return cfg_vec;
	}
	// - NQSNH::PsiUpdate will update the variational parameters of the wavefunction when given a vector of changes.
	void NQSNH::PsiUpdate(VectorXd dP)
	{
		dP = dP.array() * OptInds.array();
		for (int p = 0; p < Np; p++)
		{
			if (abs(dP(p)) < 1e-30)
			{
				dP(p) = 0;
			}
		}
		double da = dP(0);
		double dA = dP(1);
		VectorXd db = dP.segment(2, Alpha);
		VectorXd dB = dP.segment(2 + Alpha, Alpha);
		VectorXd dW = dP.tail(Alpha * Nv);
		MatrixXd dWv(Map<MatrixXd>(dW.data(), Nv, Alpha));
		dWv.transposeInPlace();
		a += da;
		A += dA;
		b += db;
		B += dB;
		Wv += dWv;
		ParamCheck();
		ParamFill();
		return;
	}
	// - NQSNH::PrepPsi will load any local configuration information used in the wavefunction constituents.
	void NQSNH::PrepPsi(Config* Cfg)
	{
		VectorXd cfg_vec = CfgRead(Cfg);
		Theta = bv + (Wm * cfg_vec);
		Nsq = cfg_vec.array().square();
		ThetaP = Theta;
		NsqP = Nsq;
		return;
	}
	// - NQSNH::PsiCfgUpdate will update any local configuration information after a configuration change.
	void NQSNH::PsiCfgUpdate() // This version uses the stored alternate local information.
	{
		Theta = ThetaP;
		Nsq = NsqP;
		return;
	}
	// - NQSNH::PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
	// -- appropriate local information for one configuration is loaded into the Ansatz.
	double NQSNH::PsiRatio(Diff diff)
	{
		double Ratio = 1;
		VectorXd ThetaShift = VectorXd::Zero(Nh); // May need to redo this initialisation.
		VectorXd NsqShift = VectorXd::Zero(Nv);
		for (int d = 0; d < diff.num; d++)
		{
			Ratio *= exp(av(diff.pos[d]) * diff.val[d]);
			ThetaShift += (Wm.col(diff.pos[d]) * diff.val[d]);
			NsqShift(diff.pos[d]) = ((2 * sqrt(Nsq(diff.pos[d])) + diff.val[d]) * diff.val[d]);
			Ratio *= exp(Av(diff.pos[d]) * NsqShift(diff.pos[d]));
		}
		NsqP = Nsq + NsqShift;
		ThetaP = Theta + ThetaShift;
		VectorXd Trace = NHTrace(Theta);
		VectorXd TraceP = NHTrace(ThetaP);
		VectorXd TraceRatio = TraceP.array() / Trace.array();
		Ratio *= TraceRatio.prod();
		if (isnan(Ratio) || isinf(Ratio))
		{
			Ratio = 0;
		}
		return Ratio;
	}
	// - NQSNH::LogDeriv will return the logarithmic derivatives of the wavefunction w.r.t. its parameters.
	VectorXd NQSNH::LogDeriv(Config* Cfg) const
	{
		VectorXd dLogp = VectorXd::Zero(Np);
		VectorXd cfg_vec = CfgRead(Cfg);
		int Ntr = g_->Ntranslate();
		if (OptInds(0))
		{
			dLogp(0) = cfg_vec.sum();
		}
		if (OptInds(1))
		{
			dLogp(1) = Nsq.sum();
		}
		VectorXd dTheta = dT_NHTrace(Theta);
		VectorXd dB = dB_NHTrace(Theta);
		for (int al = 0; al < Alpha; al++)
		{
			if (OptInds(al + 2))
			{
				VectorXd dThetaA = dTheta.segment(al * Ntr, Ntr);
				dLogp(al + 2) = dThetaA.sum();
			}
			if (OptInds(al + Alpha + 2))
			{
				VectorXd dBA = dB.segment(al * Ntr, Ntr);
				dLogp(al + 2 + Alpha) = dBA.sum();
			}
			for (int b = 0; b < Ntr; b++)
			{
				vector<int> Bond = g_->BondSearch(b);
				int TInd = b + al * Ntr;
				for (int n = 0; n < Nv; n++)
				{
					int PInd = 2 + (2 * Alpha) + n + (al * Nv);
					int VInd = Bond[n];
					if (OptInds(PInd))
					{
						dLogp(PInd) += (cfg_vec(VInd) * dTheta(TInd));
					}
				}
			}
		}
		for (int np = 0; np < Np; np++)
		{
			if (isnan(dLogp(np)) || isinf(dLogp(np)))
			{
				dLogp(np) = 0;
			}
		}
		return dLogp;
	}

	// <<<<<<------------ NQSMH functions ----------->>>>>>
	// Constructors for the multiplon-holon NQS:
		// - Random initialisation with starting values of each parameter type.
	NQSMH::NQSMH(Hilbert* hlb, Graph* grp, int HUDen, vector<double> StartParams, vector<double> noise)
	{
		NetworkInit(hlb, grp, HUDen);
		double nmag = noise[0]; // Noise vector should be populated in order: nmag, nphs.
		random_device rd;
		default_random_engine rndgen(rd());
		uniform_real_distribution<double> noisedist(-nmag, nmag);
		// Initialise random parameters, no random phase initially.
		a = StartParams[0] * (1 + noisedist(rndgen));
		A = StartParams[1] * (1 + noisedist(rndgen));
		for (int a = 0; a < Alpha; a++)
		{
			BH(a) = StartParams[2] * (1 + noisedist(rndgen));
			BM(a) = StartParams[3] * (1 + noisedist(rndgen));
			for (int n = 0; n < Nv; n++)
			{
				Wv(a, n) = StartParams[4] * (1 + noisedist(rndgen));
				Xv(a, n) = StartParams[5] * (1 + noisedist(rndgen));
			}
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
	NQSMH::NQSMH(Hilbert* hlb, Graph* grp, int HUDen, VectorXd ParameterList)
	{
		NetworkInit(hlb, grp, HUDen);
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
	// - NQSMH::NetworkInit will initialise the visible and hidden units of the RBM.
	void NQSMH::NetworkInit(Hilbert* hlb, Graph* grp, int HUDen)
	{
		h_ = hlb;
		g_ = grp;
		if (h_->Type() == 'f')
		{
			cerr << "NQSMH is not compatible with fermionic systems." << endl;
			std::abort();
		}
		if (h_->SysSize() != g_->Nsite())
		{
			cerr << "Hilbert and Graph lattice size mismatch." << endl;
			std::abort();
		}
		Alpha = HUDen;
		HDim = h_->SiteDim();
		Nv = h_->SysSize();
		Nh = Alpha * g_->Ntranslate();
		Np = 2 + 2 * Alpha * (1 + Nv);
		av.resize(Nv);
		Av.resize(Nv);
		BH.resize(Alpha);
		BM.resize(Alpha);
		BHv.resize(Nh);
		BMv.resize(Nh);
		Wv.resize(Alpha, Nv);
		Xv.resize(Alpha, Nv);
		Wm.resize(Nh, Nv);
		Xm.resize(Nh, Nv);
		ThetaH.resize(Nh);
		ThetaM.resize(Nh);
		Nsq.resize(Nv);
		Hv.resize(Nv);
		Mv.resize(Nv);
		return;
	}
	// - NQSMH::ParamFill will populate the vectors and matrices given relevant parameters.
	void NQSMH::ParamFill()
	{
		int Ntr = g_->Ntranslate();
		for (int n = 0; n < Nv; n++)
		{
			av(n) = a;
			Av(n) = A;
		}
		for (int h = 0; h < Alpha; h++)
		{
			for (int t = 0; t < Ntr; t++)
			{
				int HInd = t + h * Ntr;
				BHv(HInd) = BH(h);
				BMv(HInd) = BM(h);
				vector<int> Bond = g_->BondSearch(t);
				for (int n = 0; n < Nv; n++)
				{
					if (Bond[n] >= 0)
					{
						Wm(HInd, Bond[n]) = Wv(h, n);
						Xm(HInd, Bond[n]) = Xv(h, n);
					}
				}
			}
		}
		return;
	}
	// - NQSMH::ParamCheck will test the values of all parameters and ensure that they are less than the cap.
	void NQSMH::ParamCheck()
	{
		if (abs(a) > ParamCap)
		{
			a = a * ParamCap / abs(a);
		}
		if (abs(A) > ParamCap)
		{
			A = A * ParamCap / abs(A);
		}
		if (isnan(a) || isinf(a))
		{
			a = 0;
		}
		if (isnan(A) || isinf(A))
		{
			A = 0;
		}
		for (int al = 0; al < Alpha; al++)
		{
			if (abs(BH(al)) > ParamCap)
			{
				BH(al) = BH(al) * ParamCap / abs(BH(al));
			}
			if (abs(BM(al)) > ParamCap)
			{
				BM(al) = BM(al) * ParamCap / abs(BM(al));
			}
			if (isnan(BH(al)) || isinf(BH(al)))
			{
				BH(al) = 0;
			}
			if (isnan(BM(al)) || isinf(BM(al)))
			{
				BM(al) = 0;
			}
			for (int n = 0; n < Nv; n++)
			{
				if (abs(Wv(al, n)) > ParamCap)
				{
					Wv(al, n) = Wv(al, n) * ParamCap / abs(Wv(al, n));
				}
				if (isnan(Wv(al, n)) || isinf(Wv(al, n)))
				{
					Wv(al, n) = 0;
				}
				if (abs(Xv(al, n)) > ParamCap)
				{
					Xv(al, n) = Xv(al, n) * ParamCap / abs(Xv(al, n));
				}
				if (isnan(Xv(al, n)) || isinf(Xv(al, n)))
				{
					Xv(al, n) = 0;
				}
			}
		}
		return;
	}
	// Observer functions:
	// - NQSMH::Nvisible returns the number of visible sites.
	int NQSMH::Nvisible() const
	{
		return Nv;
	}
	// - NQSMH::Nhidden returns the number of hidden sites.
	int NQSMH::Nhidden() const
	{
		return Nh;
	}
	// - NQSMH::Nparam returns the number of parameters.
	int NQSMH::Nparam() const
	{
		return Np;
	}
	// - NQSMH::VarFlag will return 0 or 1 according to whether variational modification is permitted.
	bool NQSMH::VarFlag() const
	{
		return VFlag;
	}
	// - NQSMH::OptIndRead returns the individual parameter flag vector.
	VectorXd NQSMH::OptIndRead() const
	{
		return OptInds;
	}
	// - NQSMH::ParamLoad will change the existing parameters to the provided parameters.
	void NQSMH::ParamLoad(VectorXd NewParams)
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
		NewParams = NewParams.array() * OptInds.array();
		double da = NewParams(0);
		double dA = NewParams(1);
		VectorXd dBH = NewParams.segment(2, Alpha);
		VectorXd dBM = NewParams.segment(2 + Alpha, Alpha);
		VectorXd dW = NewParams.segment(2 + 2 * Alpha, Alpha * Nv);
		VectorXd dX = NewParams.tail(Alpha * Nv);
		MatrixXd dWv(Map<MatrixXd>(dW.data(), Nv, Alpha));
		MatrixXd dXv(Map<MatrixXd>(dX.data(), Nv, Alpha));
		dWv.transposeInPlace();
		dXv.transposeInPlace();
		a = da;
		A = dA;
		BH = dBH;
		BM = dBM;
		Wv = dWv;
		Xv = dXv;
		double MaxParam = max(abs(NewParams.maxCoeff()), abs(NewParams.minCoeff()));
		ParamCap = max(ParamCap, MaxParam);
		ParamCheck();
		ParamFill();
		return;
	}
	// - NQSMH::ParamList returns the parameters in the Modifier
	// -- List order is <a, A, BH, BM, Wv, Xv>.
	VectorXd NQSMH::ParamList() const
	{
		VectorXd Params;
		Params.resize(Np);
		Params(0) = a;
		Params(1) = A;
		Params.segment(2, Alpha) = BH;
		Params.segment(2 + Alpha, Alpha) = BM;
		MatrixXd Wt = Wv.transpose();
		MatrixXd Xt = Xv.transpose();
		VectorXd Wvec(Map<VectorXd>(Wt.data(), Wv.cols() * Wv.rows()));
		VectorXd Xvec(Map<VectorXd>(Xt.data(), Xv.cols() * Xv.rows()));
		Params.segment(2 + 2 * Alpha, Alpha * Nv) = Wvec;
		Params.tail(Alpha * Nv) = Wvec;
		return Params;
	}

	// Variational modification management functions:
	// - VarSwitch will alter VFlag if allowed.
	void NQSMH::VarSwitch()
	{
		VFlag = !VFlag;
		return;
	}
	// - OptIndLoad will load a Boolean vector of updated parameter flags.
	void NQSMH::OptIndLoad(VectorXd newinds)
	{
		if (newinds.size() != Np)
		{
			cerr << "New optimisation index vector does not have the correct number of entries." << endl;
			std::abort();
		}
		OptInds = newinds;
		return;
	}
	// - RndBtcSelect will randomly alter the individual parameter flags.
	// -- Default zero argument version will select 1/e of the parameters.
	void NQSMH::RndBtcSelect()
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
	void NQSMH::RndBtcSelect(double pfrac)
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
	// - NQSMH::SetParamCap will change the maximum allowed amplitude of a single parameter.
	void NQSMH::SetParamCap(double newcap)
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
	VectorXd NQSMH::MHTrace(VectorXd thetah_in, VectorXd thetam_in) const // Calculate trace over hidden units for given ThetaH/M.
	{
		VectorXd Hh = VectorXd::Zero(HDim);
		Hh(0) = 1.0;
		VectorXd Mh = VectorXd::LinSpaced(HDim, -1, (int)(HDim - 2));
		Mh(0) = 0;
		MatrixXd F = (thetah_in * Hh.transpose()) + (thetam_in * Mh.transpose());
		F = F.array().exp();
		VectorXd Tr = VectorXd::Zero(F.rows()); //  F.rowwise().sum() / HDim;
		for (int c = 0; c < F.cols(); c++)
		{
			Tr += (F.col(c) / (double)HDim);
		}
		return Tr;
	}
	VectorXd NQSMH::dTH_MHTrace(VectorXd thetah_in, VectorXd thetam_in) const// Calculate derivative w.r.t ThetaH.
	{
		VectorXd Hh = VectorXd::Zero(HDim);
		Hh(0) = 1.0;
		VectorXd Mh = VectorXd::LinSpaced(HDim, -1, (int)(HDim - 2));
		Mh(0) = 0;
		MatrixXd F = (thetah_in * Hh.transpose()) + (thetam_in * Mh.transpose());
		F = F.array().exp();
		VectorXd dF = F.col(0);
		VectorXd Tr = VectorXd::Zero(F.rows());
		for (int c = 0; c < F.cols(); c++)
		{
			Tr += F.col(c);
		}
		VectorXd dTr = dF.array() / Tr.array();
		return dTr;
	}
	VectorXd NQSMH::dTM_MHTrace(VectorXd thetah_in, VectorXd thetam_in) const// Calculate derivative w.r.t ThetaM.
	{
		VectorXd Hh = VectorXd::Zero(HDim);
		Hh(0) = 1.0;
		VectorXd Mh = VectorXd::LinSpaced(HDim, -1, (int)(HDim - 2));
		Mh(0) = 0;
		MatrixXd Mhm = VectorXd::Ones(thetah_in.size()) * Mh.transpose();
		MatrixXd F = (thetah_in * Hh.transpose()) + (thetam_in * Mh.transpose());
		F = F.array().exp();
		MatrixXd dF = F.array() * Mhm.array();
		VectorXd dFv = VectorXd::Zero(F.rows());
		VectorXd Tr = VectorXd::Zero(F.rows());
		for (int c = 0; c < F.cols(); c++)
		{
			dFv += dF.col(c);
			Tr += F.col(c);
		}
		VectorXd dTr = dFv.array() / Tr.array();
		return dTr;
	}
	// - NQSMH::CfgRead will allow for additional modification of FullCfg methods in Config.
	VectorXd NQSMH::CfgRead(Config* Cfg) const
	{
		vector<int> cfg_init = Cfg->FullCfg();
		vector<double> cfg_double(cfg_init.begin(), cfg_init.end());
		double* cfg_p = &cfg_double[0];
		Map<VectorXd> cfg_vec(cfg_p, cfg_double.size());
		return cfg_vec;
	}
	// - NQSMH::PsiUpdate will update the variational parameters of the wavefunction when given a vector of changes.
	void NQSMH::PsiUpdate(VectorXd dP)
	{
		dP = dP.array() * OptInds.array();
		for (int p = 0; p < Np; p++)
		{
			if (abs(dP(p)) < 1e-30)
			{
				dP(p) = 0;
			}
		}
		double da = dP(0);
		double dA = dP(1);
		VectorXd dBH = dP.segment(2, Alpha);
		VectorXd dBM = dP.segment(2 + Alpha, Alpha);
		VectorXd dW = dP.segment(2 + 2 * Alpha, Alpha * Nv);
		VectorXd dX = dP.tail(Alpha * Nv);
		MatrixXd dWv(Map<MatrixXd>(dW.data(), Nv, Alpha));
		MatrixXd dXv(Map<MatrixXd>(dX.data(), Nv, Alpha));
		dWv.transposeInPlace();
		dXv.transposeInPlace();
		a += da;
		A += dA;
		BH += dBH;
		BM += dBM;
		Wv += dWv;
		Xv += dXv;
		ParamCheck();
		ParamFill();
		return;
	}
	// - NQSMH::PrepPsi will load any local configuration information used in the wavefunction constituents.
	void NQSMH::PrepPsi(Config* Cfg)
	{
		VectorXd cfg_vec = CfgRead(Cfg);
		Nsq = cfg_vec.array().square();
		for (int n = 0; n < cfg_vec.size(); n++)
		{
			if (cfg_vec(n) == 0)
			{
				Hv(n) = 1;
			}
			else
			{
				Mv(n) = cfg_vec(n) - 1;
			}
		}
		ThetaH = BHv + (Wm * Hv) + (Xm * Mv);
		ThetaM = BMv + (Wm * Mv) + (Xm * Hv);
		ThetaHP = ThetaH;
		ThetaMP = ThetaM;
		NsqP = Nsq;
		HvP = Hv;
		MvP = Mv;
		return;
	}
	// - NQSMH::PsiCfgUpdate will update any local configuration information after a configuration change.
	void NQSMH::PsiCfgUpdate() // This version uses the stored alternate local information.
	{
		ThetaH = ThetaHP;
		ThetaM = ThetaMP;
		Nsq = NsqP;
		Hv = HvP;
		Mv = MvP;
		return;
	}
	// - NQSMH::PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
	// -- appropriate local information for one configuration is loaded into the Ansatz.
	double NQSMH::PsiRatio(Diff diff)
	{
		double Ratio = 1;
		VectorXd ThetaHShift = VectorXd::Zero(Nh);
		VectorXd ThetaMShift = VectorXd::Zero(Nh);
		VectorXd NsqShift = VectorXd::Zero(Nv);
		VectorXd HvShift = VectorXd::Zero(Nv);
		VectorXd MvShift = VectorXd::Zero(Nv);
		for (int d = 0; d < diff.num; d++)
		{
			Ratio *= exp(av(diff.pos[d]) * diff.val[d]);
			HvShift(diff.pos[d]) = ((Mv(diff.pos[d]) + diff.val[d]) < 0) - ((diff.val[d] < 0) * Hv(diff.pos[d]));
			MvShift(diff.pos[d]) = diff.val[d] + HvShift(diff.pos[d]);
			NsqShift(diff.pos[d]) = ((2 * sqrt(Nsq(diff.pos[d])) + diff.val[d]) * diff.val[d]);
			ThetaHShift += (HvShift(diff.pos[d]) * Wm.col(diff.pos[d]) + (MvShift(diff.pos[d]) * Xm.col(diff.pos[d])));
			ThetaMShift += (MvShift(diff.pos[d]) * Wm.col(diff.pos[d]) + (HvShift(diff.pos[d]) * Xm.col(diff.pos[d])));
			Ratio *= exp(Av(diff.pos[d]) * NsqShift(diff.pos[d]));
		}
		NsqP = Nsq + NsqShift;
		ThetaHP = ThetaH + ThetaHShift;
		ThetaMP = ThetaM + ThetaMShift;
		HvP = Hv + HvShift;
		MvP = Mv + MvShift;
		VectorXd Trace = MHTrace(ThetaH, ThetaM);
		VectorXd TraceP = MHTrace(ThetaHP, ThetaMP);
		VectorXd TraceRatio = TraceP.array() / Trace.array();
		Ratio *= TraceRatio.prod();
		if (isnan(Ratio) || isinf(Ratio))
		{
			Ratio = 0;
		}
		return Ratio;
	}
	// - NQSMH::LogDeriv will return the logarithmic derivatives of the wavefunction w.r.t. its parameters.
	VectorXd NQSMH::LogDeriv(Config* Cfg) const
	{
		VectorXd dLogp = VectorXd::Zero(Np);
		VectorXd cfg_vec = CfgRead(Cfg);
		int Ntr = g_->Ntranslate();
		if (OptInds(0))
		{
			dLogp(0) = cfg_vec.sum();
		}
		if (OptInds(1))
		{
			dLogp(1) = Nsq.sum();
		}
		VectorXd dThetaH = dTH_MHTrace(ThetaH, ThetaM);
		VectorXd dThetaM = dTM_MHTrace(ThetaH, ThetaM);
		for (int al = 0; al < Alpha; al++)
		{
			if (OptInds(al + 2))
			{
				VectorXd dThetaHA = dThetaH.segment(al * Ntr, Ntr);
				dLogp(al + 2) = dThetaHA.sum();
			}
			if (OptInds(al + 2 + Alpha))
			{
				VectorXd dThetaMA = dThetaM.segment(al * Ntr, Ntr);
				dLogp(al + 2 + Alpha) = dThetaMA.sum();
			}
			for (int b = 0; b < Ntr; b++)
			{
				vector<int> Bond = g_->BondSearch(b);
				int TInd = b + al * Ntr;
				for (int n = 0; n < Nv; n++)
				{
					int PIndW = 2 + (2 * Alpha) + n + (al * Nv);
					int PIndX = 2 + (2 * Alpha) + n + ((al + Alpha) * Nv);
					int VInd = Bond[n];
					if (OptInds(PIndW))
					{
						dLogp(PIndW) += (Hv(VInd) * dThetaH(TInd)) + (Mv(VInd) * dThetaM(TInd));
					}
					if (OptInds(PIndX))
					{
						dLogp(PIndX) += (Mv(VInd) * dThetaH(TInd)) + (Hv(VInd) * dThetaM(TInd));
					}
				}
			}
		}
		for (int np = 0; np < Np; np++)
		{
			if (isnan(dLogp(np)) || isinf(dLogp(np)))
			{
				dLogp(np) = 0;
			}
		}
		return dLogp;
	}

	// Constructors for the spin-hidden one-hot NQS:
	// - Random initialisation with starting values of each parameter type.
	NQSOH::NQSOH(Hilbert* hlb, Graph* grp, int HUDen, vector<double> StartParams, vector<double> noise)
	{
		NetworkInit(hlb, grp, HUDen);
		double nmag = noise[0]; // Noise vector should be populated in order: nmag, nphs.
		random_device rd;
		default_random_engine rndgen(rd());
		uniform_real_distribution<double> noisedist(-nmag, nmag);
		// Initialise random parameters, no random phase initially.
		for (int v = 0; v < VDim; v++)
		{
			a(v) = StartParams[0] * (1 + noisedist(rndgen));
		}
		for (int al = 0; al < Alpha; al++)
		{
			b(al) = StartParams[1] * (1 + noisedist(rndgen));
			for (int n = 0; n < (VDim * Nv); n++)
			{
				Wv(al, n) = StartParams[2] * (1 + noisedist(rndgen));
			}
		}
		// Automatically populate OptInds with ones.
		OptInds.resize(Np);
		OptInds.setOnes();
		// Starting parameter cap of 5.
		ParamCap = 5;
		ParamFill();
		// Variational flag automatically set to 1.
		VFlag = 1;

		return;
	}
	// - Initialisation with pre-existing parameters.
	NQSOH::NQSOH(Hilbert* hlb, Graph* grp, int HUDen, VectorXd ParameterList)
	{
		NetworkInit(hlb, grp, HUDen);
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
	// - NQSOH::NetworkInit will initialise the visible and hidden units of the RBM.
	void NQSOH::NetworkInit(Hilbert* hlb, Graph* grp, int HUDen)
	{
		h_ = hlb;
		g_ = grp;
		if (h_->Type() == 'f')
		{
			cerr << "NQSOH is not compatible with fermionic systems." << endl;
			std::abort();
		}
		if (h_->SysSize() != g_->Nsite())
		{
			cerr << "Hilbert and Graph lattice size mismatch." << endl;
			std::abort();
		}
		Alpha = HUDen;
		VDim = h_->SiteDim();
		if (h_->Type() == 's')
		{
			VList = VectorXd::LinSpaced(VDim, (1 - VDim) / 2, (VDim - 1) / 2).array() * (2 - (VDim % 2));
		}
		else if (h_->Type() == 'b')
		{
			VList = VectorXd::LinSpaced(VDim, 0, (VDim - 1));
		}
		Nv = h_->SysSize();
		Nh = Alpha * g_->Ntranslate();
		Np = VDim + Alpha + Alpha * VDim * Nv;
		a.resize(VDim);
		av.resize(VDim * Nv);
		b.resize(Alpha);
		bv.resize(Nh);
		Wv.resize(Alpha, VDim * Nv);
		Wm.resize(Nh, VDim * Nv);
		Theta.resize(Nh);
		OHVec.resize(VDim * Nv);
		return;
	}
	// - NQSOH::ParamFill will populate the vectors and matrices given relevant parameters.
	void NQSOH::ParamFill()
	{
		int HInd;
		int VInd;
		int Ntr = g_->Ntranslate();
		for (int v = 0; v < VDim; v++)
		{
			for (int n = 0; n < Nv; n++)
			{
				VInd = v + n * VDim;
				av(VInd) = a(v);
			}
		}
		for (int h = 0; h < Alpha; h++)
		{
			for (int t = 0; t < Ntr; t++)
			{
				HInd = t + h * Ntr;
				bv(HInd) = b(h);
				vector<int> Bond = g_->BondSearch(t);
				for (int n = 0; n < Nv; n++)
				{
					if (Bond[n] >= 0)
					{
						for (int v = 0; v < VDim; v++)
						{
							VInd = v + Bond[n] * VDim;
							Wm(HInd, VInd) = Wv(h, v + n * VDim);
						}
					}
				}
			}
		}
		return;
	}
	// - NQSOH::ParamCheck will test the values of all parameters and ensure that they are less than the cap.
	void NQSOH::ParamCheck()
	{
		for (int v = 0; v < VDim; v++)
		{
			if (abs(a(v)) > ParamCap)
			{
				a(v) = a(v) * ParamCap / abs(a(v));
			}
			if (isnan(a(v)) || isinf(a(v)))
			{
				a(v) = 0;
			}
		}
		for (int al = 0; al < Alpha; al++)
		{
			if (abs(b(al)) > ParamCap)
			{
				b(al) = b(al) * ParamCap / abs(b(al));
			}
			if (isnan(b(al)) || isinf(b(al)))
			{
				b(al) = 0;
			}
			for (int n = 0; n < (VDim * Nv); n++)
			{
				if (abs(Wv(al, n)) > ParamCap)
				{
					Wv(al, n) = Wv(al, n) * ParamCap / abs(Wv(al, n));
				}
				if (isnan(Wv(al, n)) || isinf(Wv(al, n)))
				{
					Wv(al, n) = 0;
				}
			}
		}
		return;
	}
	// Observer functions:
	// - NQSOH::Nvisible returns the number of visible sites.
	int NQSOH::Nvisible() const
	{
		return Nv;
	}
	// - NQSOH::Nhidden returns the number of hidden sites.
	int NQSOH::Nhidden() const
	{
		return Nh;
	}
	// - NQSOH::Nparam returns the number of parameters.
	int NQSOH::Nparam() const
	{
		return Np;
	}
	// - NQSOH::VarFlag will return 0 or 1 according to whether variational modification is permitted.
	bool NQSOH::VarFlag() const
	{
		return VFlag;
	}
	// - NQSOH::OptIndRead returns the individual parameter flag vector.
	VectorXd NQSOH::OptIndRead() const
	{
		return OptInds;
	}
	// - NQSOH::ParamLoad will change the existing parameters to the provided parameters.
	void NQSOH::ParamLoad(VectorXd NewParams)
	{
		if (NewParams.size() != Np)
		{
			cerr << "Number of provided parameters does not match Modifier parameter number." << endl;
			std::abort();
		}
		for (int p = 0; p < Np; p++)
		{
			if (abs(NewParams(p)) < 1e-30)
			{
				NewParams(p) = 0;
			}
		}
		NewParams = NewParams.array() * OptInds.array();
		VectorXd da = NewParams.segment(0, VDim);
		VectorXd db = NewParams.segment(VDim, Alpha);
		VectorXd dW = NewParams.tail(Alpha * VDim * Nv);
		MatrixXd dWv(Map<MatrixXd>(dW.data(), (VDim * Nv), Alpha));
		dWv.transposeInPlace();
		a = da;
		b = db;
		Wv = dWv;
		double MaxParam = max(abs(NewParams.maxCoeff()), abs(NewParams.minCoeff()));
		ParamCap = max(ParamCap, MaxParam);
		ParamCheck();
		ParamFill();
		return;
	}

	// - NQSOH::ParamList returns the parameters in the Modifier
	// -- List order is <a, b, Wv>.
	VectorXd NQSOH::ParamList() const
	{
		VectorXd Params;
		Params.resize(Np);
		Params.segment(0, VDim) = a;
		Params.segment(VDim, Alpha) = b;
		MatrixXd Wt = Wv.transpose();
		VectorXd Wvec(Map<VectorXd>(Wt.data(), Wv.cols() * Wv.rows()));
		Params.tail(Alpha * VDim * Nv) = Wvec;
		return Params;
	}

	// Variational modification management functions:
	// - NQSOH::VarSwitch will alter VFlag if allowed.
	void NQSOH::VarSwitch()
	{
		VFlag = !VFlag;
		return;
	}
	// - NQSOH::OptIndLoad will load a Boolean vector of updated parameter flags.
	void NQSOH::OptIndLoad(VectorXd newinds)
	{
		if (newinds.size() != Np)
		{
			cerr << "New optimisation index vector does not have the correct number of entries." << endl;
			std::abort();
		}
		OptInds = newinds;
		return;
	}
	// - NQSOH::RndBtcSelect will randomly alter the individual parameter flags.
	// -- Default zero argument version will select 1/e of the parameters.
	void NQSOH::RndBtcSelect()
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
	void NQSOH::RndBtcSelect(double pfrac)
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
	// - NQSOH::SetParamCap will change the maximum allowed amplitude of a single parameter.
	void NQSOH::SetParamCap(double newcap)
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

	// - NQSOH::CfgRead will allow for additional modification of FullCfg methods in Config.
	VectorXd NQSOH::CfgRead(Config* Cfg) const // Reads in configuration and outputs one-hot version.
	{
		vector<int> cfg_init = Cfg->FullCfg();
		VectorXd cfg_vec = VectorXd::Zero(cfg_init.size() * VDim);
		for (int n = 0; n < cfg_init.size(); n++)
		{
			for (int v = 0; v < VDim; v++)
			{
				if (cfg_init[n] == VList(v))
				{
					cfg_vec(VDim * n + v) = 1;
				}
			}
		}
		return cfg_vec;
	}
	// - NQSOH::PsiUpdate will update the variational parameters of the wavefunction when given a vector of changes.
	void NQSOH::PsiUpdate(VectorXd dP)
	{
		dP = dP.array() * OptInds.array();
		VectorXd da = dP.segment(0, VDim);
		VectorXd db = dP.segment(VDim, Alpha);
		VectorXd dW = dP.tail(Alpha * VDim * Nv);
		MatrixXd dWv(Map<MatrixXd>(dW.data(), (VDim * Nv), Alpha));
		dWv.transposeInPlace();
		a = da;
		b = db;
		Wv = dWv;
		ParamCheck();
		ParamFill();
		return;
	}
	// - NQSOH::PrepPsi will load any local configuration information used in the wavefunction constituents.
	void NQSOH::PrepPsi(Config* Cfg)
	{
		VectorXd cfg_vec = CfgRead(Cfg);
		Theta = bv + (Wm * cfg_vec);
		OHVec = cfg_vec;
		ThetaP = Theta;
		OHVecP = OHVec;
		return;
	}
	// - NQSOH::PsiCfgUpdate will update any local configuration information after a configuration change.
	void NQSOH::PsiCfgUpdate() // This version uses the stored alternate local information.
	{
		Theta = ThetaP;
		OHVec = OHVecP;
		return;
	}
	// - NQSOH::PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
	// -- appropriate local information for one configuration is loaded into the Ansatz.
	double NQSOH::PsiRatio(Diff diff)
	{
		double dV = VList(1) - VList(0);
		double Ratio = 1;
		int Ind0;
		int IndP;
		int SegStart;
		VectorXd IndAncilla = VectorXd::LinSpaced(VDim, 0, (int)(VDim - 1));
		VectorXd OHVecShift = VectorXd::Zero((int)(VDim * Nv));
		VectorXd OHSeg;
		for (int d = 0; d < diff.num; d++)
		{
			SegStart = (VDim * diff.pos[d]);
			OHSeg = OHVec.segment(SegStart, VDim).array() * IndAncilla.array();
			Ind0 = (int)OHSeg.sum();
			IndP = ((int)(Ind0 + (diff.val[d] / dV)) % VDim) + SegStart;
			OHVecShift(Ind0 + SegStart) = -OHVec(Ind0 + SegStart);
			OHVecShift(IndP) = 1;
		}
		VectorXd ThetaShift = Wm * OHVecShift;
		OHVecP = OHVec + OHVecShift;
		ThetaP = Theta + ThetaShift;
		VectorXd Trace = Theta.array().cosh();
		VectorXd TraceP = ThetaP.array().cosh();
		VectorXd TraceRatio = TraceP.array() / Trace.array();
		VectorXd dVA = av.array() * OHVecShift.array();
		Ratio *= exp(dVA.sum()) * TraceRatio.prod();
		if (isnan(Ratio) || isinf(Ratio))
		{
			Ratio = 0;
		}
		return Ratio;
	}
	// - NQSOH::LogDeriv will return the logarithmic derivatives of the wavefunction w.r.t. its parameters.
	VectorXd NQSOH::LogDeriv(Config* Cfg) const
	{
		VectorXd dLogp = VectorXd::Zero(Np);
		VectorXd cfg_vec = CfgRead(Cfg);
		int Ntr = g_->Ntranslate();
		for (int v = 0; v < VDim; v++)
		{
			if (OptInds(v))
			{
				for (int n = 0; n < Nv; n++)
				{
					dLogp(v) += cfg_vec(n * VDim + v);
				}
			}
		}
		VectorXd dTheta = Theta.array().tanh();
		int PSegStart;
		int CfgSegStart;
		vector<int> Bond;
		int TInd;
		VectorXd dThetaA;
		for (int al = 0; al < Alpha; al++)
		{
			if (OptInds(al + VDim))
			{
				dThetaA = dTheta.segment(al * Ntr, Ntr);
				dLogp(al + VDim) = dThetaA.sum();
			}
			for (int n = 0; n < Nv; n++)
			{
				PSegStart = VDim + Alpha + (VDim * (al * Nv + n));
				for (int v = 0; v < VDim; v++)
				{
					if (OptInds(PSegStart + v))
					{
						for (int b = 0; b < Ntr; b++)
						{
							Bond = g_->BondSearch(b);
							CfgSegStart = VDim * Bond[n];
							TInd = b + al * Ntr;
							dLogp(PSegStart + v) += (cfg_vec(CfgSegStart + v) * dTheta(TInd));
						}
					}
				}
			}
		}
		for (int np = 0; np < Np; np++)
		{
			if (isnan(dLogp(np)) || isinf(dLogp(np)))
			{
				dLogp(np) = 0;
			}
		}
		return dLogp;
	}
}

#endif