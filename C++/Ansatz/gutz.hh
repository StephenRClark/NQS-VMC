// ------------------------------------------------------------------------------------------------------------
// gutz.hh - contains the declaration of the Gutzwiller class, which is a subclass of Modifier. This variant 
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
		Gutz(Hilbert* hlb, double gfac, bool VarFlag);
		// Internal parameter organisation functions.
		// - ParamCheck will test the values of all parameters and ensure that they are less than the cap.
		void ParamCheck();
		// Observer functions:
		// - Nsite returns the number of visible sites.
		int Nsite() const;
		// - Nparam returns the number of parameters.
		int Nparam() const;
		// - VarFlag will return 0 or 1 according to whether variational modification is permitted.
		bool VarFlag() const;
		// - OptIndRead returns the individual parameter flag vector.
		VectorXd OptIndRead() const;
		// - ParamLoad will change the existing parameters to the provided parameters.
		void ParamLoad(VectorXd g_new);
		// - ParamList returns the parameters in the Modifier
		// -- List order is <a, A, b, B, Wv>.
		VectorXd ParamList() const;

		// Variational modification management functions:
		// - VarSwitch will alter VFlag if allowed.
		void VarSwitch();
		// - OptIndLoad will load a Boolean vector of updated parameter flags.
		void OptIndLoad(VectorXd newinds);
		// - RndBtcSelect will randomly alter the individual parameter flags.
		// -- Default zero argument version will select 1/e of the parameters.
		void RndBtcSelect();
		// -- Second version with fraction will select input fraction if valid.
		void RndBtcSelect(double pfrac);
		// - SetParamCap will change the maximum allowed amplitude of a single parameter.
		void SetParamCap(double newcap);

		// Wavefunction calculation functions:
		// - CfgRead will allow for additional modification of FullCfg methods in Config.
		VectorXd CfgRead(Config* Cfg) const;
		// - PsiUpdate will update the variational parameters of the wavefunction when given a vector of changes.
		void PsiUpdate(VectorXd dP);
		// - PrepPsi will load any local configuration information used in the wavefunction constituents.
		void PrepPsi(Config* Cfg);
		// - PsiCfgUpdate will update any local configuration information after a configuration change.
		void PsiCfgUpdate(); // This version uses the stored alternate local information.
		// - PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
		// -- appropriate local information for one configuration is loaded into the Ansatz.
		double PsiRatio(Diff diff);
		// - LogDeriv will return the logarithmic derivatives of the wavefunction w.r.t. its parameters.
		VectorXd LogDeriv(Config* Cfg) const;
	};
}

#endif