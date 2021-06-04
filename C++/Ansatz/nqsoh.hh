// ------------------------------------------------------------------------------------------------------------
// nqsoh.hh - contains the definitions of the NQS class with binary hidden units and one-hot visible unit 
//			  encoding which is a subclass of Modifier. This variant will assume translation invariance, and is
//			  intended for use with bosonic or spin > 1/2 systems.
// ------------------------------------------------------------------------------------------------------------

#ifndef NQS_OH_HH
#define NQS_OH_HH

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
	class NQSOH : public Modifier
	{
	private: // Properties associated with the spin-like hidden unit NQS.
		Graph* g_; // Pointer to a Graph object for translation invariance.
		Hilbert* h_; // Pointer to a Hilbert object.
		// Network parameters:
		int Nv; // Number of visible sites, scalar.
		int Nh; // Number of hidden units, scalar.
		int Alpha; // Hidden unit density / unique stack number, scalar.
		int VDim; // Visible unit dimension, scalar.
		// Variational parameters:
		VectorXd a; // Visible bias parameters, VDim x 1 vector.
		VectorXd b; // Hidden bias parameters, Alpha x 1 vector.
		MatrixXd Wv; // Hidden-visible coupling parameters, Alpha x (VDim*Nv) matrix.
		// Full networked parameter lists:
		VectorXd av; // Linear visible bias, (VDim*Nv) x 1 vector.
		VectorXd bv; // Linear hidden bias, Nh x 1 vector.	
		MatrixXd Wm; // Hidden-visible couplings, Nh x (VDim*Nv) matrix.
		// Local information:
		VectorXd VList; // Ordered list of visible values, VDim x 1 vector.
		VectorXd OHVec; // Local record of configuration in one-hot encoding, (VDim*Nv) x 1 vector.
		VectorXd Theta; // Effective angles, Nh x 1 vector.
		VectorXd OHVecP; // New proposed one-hot configuration.
		VectorXd ThetaP; // New proposed effective angles.
		// Variational management:
		double ParamCap; // Maximum parameter magnitude, scalar.
		int Np; // Number of parameters.
		bool VFlag; // Variational modification activation flag.
		VectorXd OptInds; // Individual parameter flags, Np x 1 Boolean vector.
	public:
		// Constructors for the spin-hidden one-hot NQS:
		// - Random initialisation with starting values of each parameter type.
		NQSOH(Hilbert* hlb, Graph* grp, int HUDen, vector<double> StartParams, vector<double> noise)
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
				for (int n = 0; n < (VDim*Nv); n++)
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
		NQSOH(Hilbert* hlb, Graph* grp, int HUDen, VectorXd ParameterList)
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
				VList = VectorXd::LinSpaced(VDim, (1 - VDim) / 2, (VDim - 1) / 2).array() * (2-(VDim % 2));
			}
			else if (h_->Type() == 'b')
			{
				VList = VectorXd::LinSpaced(VDim, 0, (VDim - 1));
			}
			Nv = h_->SysSize();
			Nh = Alpha * g_->Ntranslate();
			Np = VDim + Alpha + Alpha * VDim * Nv;
			a.resize(VDim);
			av.resize(VDim*Nv);
			b.resize(Alpha);
			bv.resize(Nh);
			Wv.resize(Alpha, VDim*Nv);
			Wm.resize(Nh, VDim*Nv);
			Theta.resize(Nh);
			OHVec.resize(VDim*Nv);
			return;
		}
		// - ParamFill will populate the vectors and matrices given relevant parameters.
		void ParamFill()
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
								Wm(HInd, VInd) = Wv(h, v+n*VDim);
							}							
						}
					}
				}
			}
			return;
		}
		// - ParamCheck will test the values of all parameters and ensure that they are less than the cap.
		void ParamCheck()
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
			VectorXd da = NewParams.segment(0, VDim);
			VectorXd db = NewParams.segment(VDim, Alpha);
			VectorXd dW = NewParams.tail(Alpha*VDim*Nv);
			MatrixXd dWv(Map<MatrixXd>(dW.data(), (VDim*Nv), Alpha));
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

		// - ParamList returns the parameters in the Modifier
		// -- List order is <a, b, Wv>.
		VectorXd ParamList() const
		{
			VectorXd Params;
			Params.resize(Np);
			Params.segment(0, VDim) = a;
			Params.segment(VDim, Alpha) = b;
			MatrixXd Wt = Wv.transpose();
			VectorXd Wvec(Map<VectorXd>(Wt.data(), Wv.cols() * Wv.rows()));
			Params.tail(Alpha*VDim*Nv) = Wvec;
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
		
		// - CfgRead will allow for additional modification of FullCfg methods in Config.
		VectorXd CfgRead(Config* Cfg) const // Reads in configuration and outputs one-hot version.
		{
			vector<int> cfg_init = Cfg->FullCfg();
			VectorXd cfg_vec = VectorXd::Zero(cfg_init.size() * VDim);
			for (int n = 0; n < cfg_init.size(); n++)
			{
				for (int v = 0; v < VDim; v++)
				{
					if (cfg_init[n] == VList(v))
					{
						cfg_vec(VDim*n+v) = 1;
					}
				}
			}
			return cfg_vec;
		}
		// - PsiUpdate will update the variational parameters of the wavefunction when given a vector of changes.
		void PsiUpdate(VectorXd dP)
		{
			dP = dP.array() * OptInds.array();
			VectorXd da = dP.segment(0, VDim);
			VectorXd db = dP.segment(VDim, Alpha);
			VectorXd dW = dP.tail(Alpha*VDim*Nv);
			MatrixXd dWv(Map<MatrixXd>(dW.data(), (VDim*Nv), Alpha));
			dWv.transposeInPlace();
			a = da;
			b = db;
			Wv = dWv;
			ParamCheck();
			ParamFill();
			return;
		}
		// - PrepPsi will load any local configuration information used in the wavefunction constituents.
		void PrepPsi(Config* Cfg)
		{
			VectorXd cfg_vec = CfgRead(Cfg);
			Theta = bv + (Wm * cfg_vec);
			OHVec = cfg_vec;
			ThetaP = Theta;
			OHVecP = OHVec;
			return;
		}
		// - PsiCfgUpdate will update any local configuration information after a configuration change.
		void PsiCfgUpdate() // This version uses the stored alternate local information.
		{
			Theta = ThetaP;
			OHVec = OHVecP;
			return;
		}
		// - PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
		// -- appropriate local information for one configuration is loaded into the Ansatz.
		double PsiRatio(Diff diff)
		{
			double dV = VList(1) - VList(0);
			double Ratio = 1;
			int Ind0;
			int IndP;
			int SegStart;
			VectorXd IndAncilla = VectorXd::LinSpaced(VDim, 0, (int)(VDim-1));
			VectorXd OHVecShift = VectorXd::Zero((int)(VDim*Nv));
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
		// - LogDeriv will return the logarithmic derivatives of the wavefunction w.r.t. its parameters.
		VectorXd LogDeriv(Config* Cfg) const
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
				if (OptInds(al+VDim))
				{
					dThetaA = dTheta.segment(al*Ntr, Ntr);
					dLogp(al+VDim) = dThetaA.sum();
				}
				for (int n = 0; n < Nv; n++)
				{
					PSegStart = VDim + Alpha + (VDim * (al * Nv + n));					
					for (int v = 0; v < VDim; v++)
					{
						if (OptInds(PSegStart+v))
						{
							for (int b = 0; b < Ntr; b++)
							{
								Bond = g_->BondSearch(b);
								CfgSegStart = VDim * Bond[n];
								TInd = b + al * Ntr;
								dLogp(PSegStart+v) += (cfg_vec(CfgSegStart+v) * dTheta(TInd));
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
	};
}

#endif