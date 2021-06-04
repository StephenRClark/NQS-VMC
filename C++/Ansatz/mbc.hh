// ------------------------------------------------------------------------------------------------------------
// mbc.hh - contains the definitions of the doublon-holon many body correlator class, which is a subclass of 
//			Modifier. This variant will assume translation invariance, and is intended for use with bosonic 
//			systems.
// ------------------------------------------------------------------------------------------------------------

#ifndef MBC_HH
#define MBC_HH

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
	class NNMB : public Modifier
	{
	private: // Properties associated with the bosonic doublon-holon many body correlator.
		Graph* g_; // Pointer to a Graph object for translation invariance.
		Hilbert* h_; // Pointer to a Hilbert object.
		// System parameters:
		int N; // Number of sites, scalar.
		int z; // Coordination number of the lattice, scalar;
		vector<vector<int>> NNList; // Lists of nearest neighbours, N x z entries.
		double Nmean; // Mean number of particles per site, scalar.
		// Variational parameters:
		double gmb; // Many body correlation factor, scalar.
		// Local information:
		VectorXd dN; // Local density fluctuations, N x 1 vector.
		VectorXd dNP; // New proposed local density fluctuations.
		VectorXd Xi; // Local many body operator values, N x 1 vector.
		VectorXd XiP; // New proposed local many body operator values.
		// Variational management:
		double ParamCap; // Maximum parameter magnitude, scalar.
		int Np = 1; // Number of parameters.
		bool VFlag; // Variational modification activation flag.
		VectorXd OptInds; // Individual parameter flags, Np x 1 Boolean vector.
	public:
		// Constructors for the bosonic many body correlator:
		// - Random initialisation with starting values of each parameter type.
		NNMB(Hilbert* hlb, Graph* grp, double G0, vector<double> noise)
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
		NNMB(Hilbert* hlb, Graph* grp, VectorXd ParameterList)
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
		// - ListInit will initialise the nearest neighbour index lists.
		void ListInit(Hilbert* hlb, Graph* grp)
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
		// - ParamCheck will test the values of all parameters and ensure that they are less than the cap.
		void ParamCheck()
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
		// - Nsite returns the number of visible sites.
		int Nsite() const
		{
			return N;
		}
		// - Nparam returns the number of parameters.
		int Nparam() const
		{
			return 1;
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
		void ParamLoad(VectorXd g_new)
		{
			gmb = g_new(0);
			ParamCheck();
			return;
		}
		// - ParamList returns the parameters in the Modifier.
		VectorXd ParamList() const
		{
			VectorXd gv = VectorXd::Zero(1);
			gv(0) = gmb;
			return gv;
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
		}
		// - RndBtcSelect will randomly alter the individual parameter flags.
		// -- Default zero argument version will select 1/e of the parameters.
		void RndBtcSelect()
		{
			return; // Only one parameter, and a minimum of one is needed.
		}
		// -- Second version with fraction will select input fraction if valid.
		void RndBtcSelect(double pfrac)
		{
			return; // Only one parameter, and a minimum of one is needed.
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
			gmb += dP(0);
			ParamCheck();
			return;
		}
		// - PrepPsi will load any local configuration information used in the wavefunction constituents.
		void PrepPsi(Config* Cfg)
		{
			VectorXd cfg_vec = CfgRead(Cfg);
			Nmean = round((cfg_vec.sum() / N));
			dN = cfg_vec.array() - Nmean;
			for (int n = 0; n < N; n++)
			{
				if ((dN(n)==1) || (dN(n)==-1))
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
		// - PsiCfgUpdate will update any local configuration information after a configuration change.
		void PsiCfgUpdate() // This version uses the stored alternate local information.
		{
			Xi = XiP; dN = dNP;
			return;
		}
		// - PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
		// -- appropriate local information for one configuration is loaded into the Ansatz.
		double PsiRatio(Diff diff)
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
		// - LogDeriv will return the logarithmic derivatives of the wavefunction w.r.t. its parameters.
		VectorXd LogDeriv(Config* Cfg) const
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
	};
}

#endif