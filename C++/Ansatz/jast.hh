// ------------------------------------------------------------------------------------------------------------
// jast.hh - contains the definitions of the Jastrow class, which is a subclass of Modifier. This variant will 
//			 assume translation invariance, and is intended for use with bosonic systems.
// ------------------------------------------------------------------------------------------------------------

#ifndef JASTROW_HH
#define JASTROW_HH

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
	class Jast : public Modifier
	{
	private: // Properties associated with the bosonic Jastrow state.
		Graph* g_; // Pointer to a Graph object for translation invariance.
		Hilbert* h_; // Pointer to a Hilbert object.
		// System parameters:
		int N; // Number of sites, scalar.
		// Variational parameters:
		VectorXd JsVar; // Jastrow variables, Np x 1 vector.
		MatrixXd Js; // Jastrow factor matrix, N x N matrix.
		MatrixXd JsI; // Jastrow parameter index matrix, N x N matrix.
		// Local information:
		VectorXd Tj; // Local Jastrow factor vector, N x 1 vector.
		VectorXd TjP; // New proposed local Jastrow factors.
		// Variational management:
		double ParamCap; // Maximum parameter magnitude, scalar.
		int Np; // Number of parameters.
		bool VFlag; // Variational modification activation flag.
		VectorXd OptInds; // Individual parameter flags, Np x 1 Boolean vector.
	public:
		// Constructors for the bosonic Jastrow state:
		// - Random initialisation with starting values of each parameter type.
		Jast(Hilbert* hlb, Graph* grp, double J0, vector<double> noise)
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
		Jast(Hilbert* hlb, Graph* grp, VectorXd ParameterList)
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
		// - MatrixInit will initialise the Jastrow parameter index matrix.
		void MatrixInit(Hilbert* hlb, Graph* grp)
		{
			h_ = hlb;
			g_ = grp;
			if ((h_->Type() == 'f')||(h_->Type() == 's'))
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
								JsI(Ind1,Ind2) = pcount;
								JsI(Ind2,Ind1) = pcount;
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
		// - ParamFill will populate the vectors and matrices given relevant parameters.
		void ParamFill()
		{
			for (int n = 0; n < N; n++)
			{
				for (int m = 0; m < N; m++)
				{
					if (JsI(n, m) >= 0)
					{
						Js(n, m) = JsVar((int)JsI(n, m));
					}
				}
			}
			return;
		}
		// - ParamCheck will test the values of all parameters and ensure that they are less than the cap.
		void ParamCheck()
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
		// - Nsite returns the number of visible sites.
		int Nsite() const
		{
			return N;
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

		// - ParamList returns the parameters in the Modifier.
		VectorXd ParamList() const
		{
			return JsVar;
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
			JsVar += dP;
			ParamCheck();
			ParamFill();
			return;
		}
		// - PrepPsi will load any local configuration information used in the wavefunction constituents.
		void PrepPsi(Config* Cfg)
		{
			VectorXd cfg_vec = CfgRead(Cfg);
			Tj = Js * cfg_vec;
			TjP = Tj;
			return;
		}
		// - PsiCfgUpdate will update any local configuration information after a configuration change.
		void PsiCfgUpdate() // This version uses the stored alternate local information.
		{
			Tj = TjP;
			return;
		}
		// - PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
		// -- appropriate local information for one configuration is loaded into the Ansatz.
		double PsiRatio(Diff diff)
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
			Ratio *= exp(-DeltaMat.sum()/2);
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
	};
}

#endif