// ------------------------------------------------------------------------------------------------------------
// nqshdv.hh - contains the definitions of the NQS class with multinomial hidden units, which is a subclass of 
//			   Modifier. This variant will assume translation invariance, and is intended for use with bosonic
//			   systems.
// ------------------------------------------------------------------------------------------------------------

#ifndef NQS_HDV_HH
#define NQS_HDV_HH

#include <vector>
#include <random>
#include <Eigen/Dense>
#include "Graph/Graph.hh"
#include "Hilbert/Hilbert.hh"
#include "Modifier.hh"

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	class NQSSH : public Modifier
	{
	private: // Properties associated with the spin-like hidden unit NQS.
		Graph* g_; // Pointer to a Graph object for translation invariance.
		Hilbert* h_; // Pointer to a Hilbert object.
		// Network parameters:
		int Nv; // Number of visible sites, scalar.
		int Nh; // Number of hidden units, scalar.
		int Alpha; // Hidden unit density / unique stack number, scalar.
		int HDim; // Hidden unit dimension, scalar.
		// Variational parameters:
		double a; // Linear visible bias parameter, scalar.
		double A; // Non-linear visible bias parameter, scalar.
		VectorXd b; // Linear hidden bias parameters, Alpha x 1 vector.
		VectorXd B; // Non-linear hidden bias parameters, Alpha x 1 vector.
		MatrixXd Wv; // Hidden-visible coupling parameters, Alpha x Nv matrix.
		// Full networked parameter lists:
		VectorXd av; // Linear visible bias, Nv x 1 vector.
		VectorXd Av; // Non-linear visible bias, Nv x 1 vector.		
		VectorXd bv; // Linear hidden bias, Nh x 1 vector.
		VectorXd Bv; // Non-linear hidden bias, Nh x 1 vector.		
		MatrixXd Wm; // Hidden-visible couplings, Nh x Nv matrix.
		// Local information:
		VectorXd Nsq; // Vector of squared visible occupancies, Nv x 1 vector.
		VectorXd Theta; // Effective angles, Nh x 1 vector.
		VectorXd NsqP; // New proposed squared visible occupancies.
		VectorXd ThetaP; // New proposed effective angles.
		// Variational management:
		double ParamCap; // Maximum parameter magnitude, scalar.
		int Np; // Number of parameters.
		bool VFlag; // Variational modification activation flag.
		VectorXd OptInds; // Individual parameter flags, Np x 1 Boolean vector.
	public:
		// Constructors for the spin-hidden NQS:
		// - Random initialisation with starting values of each parameter type.
		NQSSH(Hilbert* hlb, Graph* grp, int HUDen, vector<double> StartParams, vector<double> noise)
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
		NQSSH(Hilbert* hlb, Graph* grp, int HUDen, VectorXd ParameterList)
		{
			NetworkInit(hlb, grp, HUDen);
			// Starting parameter cap of 5, or maximum given value.
			double MaxParam = max(abs(ParameterList.maxCoeff()), abs(ParameterList.minCoeff()));
			// Variational flag automatically set to 1.
			VFlag = 1;			
			// Automatically populate OptInds with ones.
			OptInds.resize(Np);
			OptInds.setOnes();
			ParamCap = max(5.0,MaxParam);
			ParamLoad(ParameterList);
			return;
		}
		// Internal parameter organisation functions.
		// - NetworkInit will initialise the visible and hidden units of the RBM.
		void NetworkInit(Hilbert* hlb, Graph* grp, int HUDen)
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
		// - ParamFill will populate the vectors and matrices given relevant parameters.
		void ParamFill()
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
		// - ParamCheck will test the values of all parameters and ensure that they are less than the cap.
		void ParamCheck()
		{
			if (abs(a) > ParamCap)
			{
				a = a * ParamCap / abs(a);
			}
			if (abs(A) > ParamCap)
			{
				A = A * ParamCap / abs(A);
			}
			if (isnan(a)||isinf(a))
			{
				a = 0;
			}
			if (isnan(A)||isinf(A))
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
						Wv(al,n) = Wv(al,n) * ParamCap / abs(Wv(al,n));
					}
					if (isnan(Wv(al, n)) || isinf(Wv(al, n)))
					{
						Wv(al,n) = 0;
					}
				}
			}
			return;
		}
		// Observer functions:
		// - Nvisible returns the number of visible sites.
		int Nvisible() const
		{
			return Nv;
		}
		// - Nhidden returns the number of hidden sites.
		int Nhidden() const
		{
			return Nh;
		}
		// - Nparam returns the number of parameters.
		int Nparam() const
		{
			return Np;
		}
		// - VarFlag will return 0 or 1 according to whether variational modification is permitted.
		bool VarFlag() const
		{
			return VFlag;
		}
		// - OptIndRead returns the individual parameter flag vector.
		VectorXd OptIndRead() const
		{
			return OptInds;
		}
		// - ParamLoad will change the existing parameters to the provided parameters.
		void ParamLoad(VectorXd NewParams)
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
			VectorXd db = NewParams.segment(2,Alpha);
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
		VectorXd ParamList() const 
		{
			VectorXd Params;
			Params.resize(Np);
			Params(0) = a;
			Params(1) = A;
			Params.segment(2, Alpha) = b;
			Params.segment(2 + Alpha,Alpha) = B;
			MatrixXd Wt = Wv.transpose();
			VectorXd Wvec(Map<VectorXd>(Wt.data(), Wv.cols() * Wv.rows()));
			Params.tail(Alpha * Nv) = Wvec;
			return Params;
		}

		// Variational modification management functions:
		// - VarSwitch will alter VFlag if allowed.
		void VarSwitch()
		{
			VFlag = !VFlag;
			return;
		}
		// - OptIndLoad will load a Boolean vector of updated parameter flags.
		void OptIndLoad(VectorXd newinds)
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
		void RndBtcSelect()
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
		void RndBtcSelect( double pfrac)
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
		// - SetParamCap will change the maximum allowed amplitude of a single parameter.
		void SetParamCap(double newcap)
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
		VectorXd SHTrace(VectorXd theta_in) // Calculate trace over hidden units for given Theta.
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
		VectorXd dT_SHTrace(VectorXd theta_in) const// Calculate derivative w.r.t Theta.
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
		VectorXd dB_SHTrace(VectorXd theta_in) const // Calculate derivative w.r.t B.
		{
			VectorXd h = VectorXd::LinSpaced(HDim, (1 - HDim)/2, (HDim - 1)/2);
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
		// - CfgRead will allow for additional modification of FullCfg methods in Config.
		VectorXd CfgRead(Config* Cfg) const
		{
			vector<int> cfg_init = Cfg->FullCfg();
			vector<double> cfg_double(cfg_init.begin(), cfg_init.end());
			double* cfg_p = &cfg_double[0];
			Map<VectorXd> cfg_vec(cfg_p, cfg_double.size());
			return cfg_vec;
		}
		// - PsiUpdate will update the variational parameters of the wavefunction when given a vector of changes.
		void PsiUpdate(VectorXd dP)
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
		// - PrepPsi will load any local configuration information used in the wavefunction constituents.
		void PrepPsi(Config* Cfg)
		{
			VectorXd cfg_vec = CfgRead(Cfg);
			Theta = bv + (Wm * cfg_vec);
			Nsq = cfg_vec.array().square();
			ThetaP = Theta;
			NsqP = Nsq;
			return;
		}
		// - PsiCfgUpdate will update any local configuration information after a configuration change.
		void PsiCfgUpdate() // This version uses the stored alternate local information.
		{
			Theta = ThetaP;
			Nsq = NsqP;
			return;
		}
		// - PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
		// -- appropriate local information for one configuration is loaded into the Ansatz.
		double PsiRatio(Diff diff)
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
		// - LogDeriv will return the logarithmic derivatives of the wavefunction w.r.t. its parameters.
		VectorXd LogDeriv(Config* Cfg) const
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
	};

	class NQSNH : public Modifier
	{
	private: // Properties associated with the spin-like hidden unit NQS.
		Graph* g_; // Pointer to a Graph object for translation invariance.
		Hilbert* h_; // Pointer to a Hilbert object.
		// Network parameters:
		int Nv; // Number of visible sites, scalar.
		int Nh; // Number of hidden units, scalar.
		int Alpha; // Hidden unit density / unique stack number, scalar.
		int HDim; // Hidden unit dimension, scalar.
		// Variational parameters:
		double a; // Linear visible bias parameter, scalar.
		double A; // Non-linear visible bias parameter, scalar.
		VectorXd b; // Linear hidden bias parameters, Alpha x 1 vector.
		VectorXd B; // Non-linear hidden bias parameters, Alpha x 1 vector.
		MatrixXd Wv; // Hidden-visible coupling parameters, Alpha x Nv matrix.
		// Full networked parameter lists:
		VectorXd av; // Linear visible bias, Nv x 1 vector.
		VectorXd Av; // Non-linear visible bias, Nv x 1 vector.		
		VectorXd bv; // Linear hidden bias, Nh x 1 vector.
		VectorXd Bv; // Non-linear hidden bias, Nh x 1 vector.		
		MatrixXd Wm; // Hidden-visible couplings, Nh x Nv matrix.
		// Local information:
		VectorXd Nsq; // Vector of squared visible occupancies, Nv x 1 vector.
		VectorXd Theta; // Effective angles, Nh x 1 vector.
		VectorXd NsqP; // New proposed squared visible occupancies.
		VectorXd ThetaP; // New proposed effective angles.
		// Variational management:
		double ParamCap; // Maximum parameter magnitude, scalar.
		int Np; // Number of parameters.
		bool VFlag; // Variational modification activation flag.
		VectorXd OptInds; // Individual parameter flags, Np x 1 Boolean vector.
	public:
		// Constructors for the spin-hidden NQS:
		// - Random initialisation with starting values of each parameter type.
		NQSNH(Hilbert* hlb, Graph* grp, int HUDen, vector<double> StartParams, vector<double> noise)
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
		NQSNH(Hilbert* hlb, Graph* grp, int HUDen, VectorXd ParameterList)
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
		// - NetworkInit will initialise the visible and hidden units of the RBM.
		void NetworkInit(Hilbert* hlb, Graph* grp, int HUDen)
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
		// - ParamFill will populate the vectors and matrices given relevant parameters.
		void ParamFill()
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
		// - ParamCheck will test the values of all parameters and ensure that they are less than the cap.
		void ParamCheck()
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
		// - Nvisible returns the number of visible sites.
		int Nvisible() const
		{
			return Nv;
		}
		// - Nhidden returns the number of hidden sites.
		int Nhidden() const
		{
			return Nh;
		}
		// - Nparam returns the number of parameters.
		int Nparam() const
		{
			return Np;
		}
		// - VarFlag will return 0 or 1 according to whether variational modification is permitted.
		bool VarFlag() const
		{
			return VFlag;
		}
		// - OptIndRead returns the individual parameter flag vector.
		VectorXd OptIndRead() const
		{
			return OptInds;
		}
		// - ParamLoad will change the existing parameters to the provided parameters.
		void ParamLoad(VectorXd NewParams)
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
		VectorXd ParamList() const
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
		// - VarSwitch will alter VFlag if allowed.
		void VarSwitch()
		{
			VFlag = !VFlag;
			return;
		}
		// - OptIndLoad will load a Boolean vector of updated parameter flags.
		void OptIndLoad(VectorXd newinds)
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
		void RndBtcSelect()
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
		void RndBtcSelect(double pfrac)
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
		// - SetParamCap will change the maximum allowed amplitude of a single parameter.
		void SetParamCap(double newcap)
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
		VectorXd NHTrace(VectorXd theta_in) // Calculate trace over hidden units for given Theta.
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
		VectorXd dT_NHTrace(VectorXd theta_in) const// Calculate derivative w.r.t Theta.
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
		VectorXd dB_NHTrace(VectorXd theta_in) const // Calculate derivative w.r.t B.
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
		// - CfgRead will allow for additional modification of FullCfg methods in Config.
		VectorXd CfgRead(Config* Cfg) const
		{
			vector<int> cfg_init = Cfg->FullCfg();
			vector<double> cfg_double(cfg_init.begin(), cfg_init.end());
			double* cfg_p = &cfg_double[0];
			Map<VectorXd> cfg_vec(cfg_p, cfg_double.size());
			return cfg_vec;
		}
		// - PsiUpdate will update the variational parameters of the wavefunction when given a vector of changes.
		void PsiUpdate(VectorXd dP)
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
		// - PrepPsi will load any local configuration information used in the wavefunction constituents.
		void PrepPsi(Config* Cfg)
		{
			VectorXd cfg_vec = CfgRead(Cfg);
			Theta = bv + (Wm * cfg_vec);
			Nsq = cfg_vec.array().square();
			ThetaP = Theta;
			NsqP = Nsq;
			return;
		}
		// - PsiCfgUpdate will update any local configuration information after a configuration change.
		void PsiCfgUpdate() // This version uses the stored alternate local information.
		{
			Theta = ThetaP;
			Nsq = NsqP;
			return;
		}
		// - PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
		// -- appropriate local information for one configuration is loaded into the Ansatz.
		double PsiRatio(Diff diff)
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
		// - LogDeriv will return the logarithmic derivatives of the wavefunction w.r.t. its parameters.
		VectorXd LogDeriv(Config* Cfg) const
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
	};
}

#endif