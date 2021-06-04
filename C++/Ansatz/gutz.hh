// ------------------------------------------------------------------------------------------------------------
// gutz.hh - contains the definitions of the Gutzwiller class, which is a subclass of Modifier. This variant 
//			 is intended for use with bosonic systems.
// ------------------------------------------------------------------------------------------------------------

#ifndef GUTWILLER_HH
#define GUTZWILLER_HH

#include <vector>
#include <random>
#include <Eigen/Dense>
#include "Hilbert/Hilbert.hh"
#include "Modifier.hh"

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	class Gutz : public Modifier
	{
	private: // Properties associated with the spin-like hidden unit NQS.
		Hilbert* h_; // Pointer to a Hilbert object.
		// System parameters:
		int N; // Number of visible sites, scalar.		
		// Variational parameters:
		double G; // Suppression factor.
		// Local information:
		double Nmean; // Mean particle density of sites - set at initialisation.
		VectorXd dN; // Vector of local density fluctuations from mean, Nv x 1 vector.
		VectorXd dNP; // Proposed density fluctuations from mean, Nv x 1 vector.
		// Variational management:
		double ParamCap = 100; // Maximum parameter magnitude, scalar.
		int Np = 1; // Number of parameters.
		bool VFlag; // Variational modification activation flag.
		VectorXd OptInds = VectorXd::Ones(1); // Individual parameter flags, Np x 1 Boolean vector.
	public:
		// Constructors for the bosonic Gutzwiller factor:
		// - Initialise with the input Gutzwiller factor gfac.
		Gutz(Hilbert* hlb, double gfac, bool VarFlag)
		{
			if ((h_->Type() == 'f')||(h_->Type() == 's'))
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
		// - ParamCheck will test the values of all parameters and ensure that they are less than the cap.
		void ParamCheck()
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
		void ParamLoad(VectorXd g_new)
		{
			G = g_new(0);
			ParamCheck();
			return;
		}

		// - ParamList returns the parameters in the Modifier
		// -- List order is <a, A, b, B, Wv>.
		VectorXd ParamList() const
		{
			VectorXd Params(1);
			Params(0) = G;
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
			G += dP(0);
			ParamCheck();
			return;
		}
		// - PrepPsi will load any local configuration information used in the wavefunction constituents.
		void PrepPsi(Config* Cfg)
		{
			VectorXd cfg_vec = CfgRead(Cfg);
			Nmean = cfg_vec.sum() / N;
			dN = cfg_vec.array() - Nmean;
			dNP = dN;
			return;
		}
		// - PsiCfgUpdate will update any local configuration information after a configuration change.
		void PsiCfgUpdate() // This version uses the stored alternate local information.
		{
			dN = dNP;
			return;
		}
		// - PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
		// -- appropriate local information for one configuration is loaded into the Ansatz.
		double PsiRatio(Diff diff)
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
		// - LogDeriv will return the logarithmic derivatives of the wavefunction w.r.t. its parameters.
		VectorXd LogDeriv(Config* Cfg) const
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
	};
}

#endif