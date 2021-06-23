// ------------------------------------------------------------------------------------------------------------
// becr.hh - contains the definitions of the Bose Einstein Condensate Reference state, which in its current 
//			 form does not support variational modification. Intended for use with bosonic systems.
// ------------------------------------------------------------------------------------------------------------

#ifndef BECR_HH
#define BECR_HH

#include <vector>
#include <cmath>
#include <random>
#include <Eigen/Dense>
#include "Graph/Graph.hh"
#include "Hilbert/Hilbert.hh"
#include "Reference.hh"

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	class BECR : public Reference
	{
	private:
		Hilbert* h_; // Pointer to a Hilbert object.
		int N; // Number of sites.
		int Nb; // Number of bosons.
		int NbP; // Proposed number of bosons.
		VectorXd Occ; // Vector of on-site occupations, N x 1 vector.
		VectorXd SPO; // Single particle orbital amplitudes, N x 1 vector.
		VectorXd OccP; // Vector of proposed on-site occupations, N x 1 vector.
	public:
		// Constructors for the Bose-Einstein Condensate Reference:
		// - Initialisation with a single particle Hamiltonian matrix.
		BECR(Hilbert* hlb, MatrixXd SPH);
		// - Initialisation with a single particle orbital vector.
		BECR(Hilbert* hlb, VectorXd Orbital);
		// Observer functions for extracting information: (Variation disabled - all return zero).
		// - Nparam will return the number of variational parameters.
		int Nparam() const;
		// - VarFlag will return 0 or 1 according to whether variational modification is permitted.
		bool VarFlag() const;
		// - OptIndRead will return individual parameter flags for variational modification.
		VectorXd OptIndRead() const;
		// - ParamList will return a vector of the variational parameters in the Reference.
		VectorXd ParamList() const;

		// Variational modification management functions: (Variation disabled - all return message)
		// - VarSwitch will alter VFlag if allowed.
		void VarSwitch();
		// - OptIndLoad will load a vector of updated parameter flags (values 0/1).
		void OptIndLoad(VectorXd newinds);
		// - RndBtcSelect will randomly alter the individual parameter flags.
		// -- Default zero argument version will select 1/e of the parameters.
		void RndBtcSelect();
		// -- Second version with fraction will select input fraction if valid.
		void RndBtcSelect(double pfrac);
		// - SetParamCap will change the maximum allowed amplitude of a single parameter.
		void SetParamCap(double newcap);
		// - ParamLoad will change the existing parameters to the provided parameters.
		void ParamLoad(VectorXd NewParams);

		// Wavefunction calculation functions:


		// - CfgRead will allow for additional modification of FullCfg methods in Config.
		VectorXd CfgRead(Config* Cfg) const;
		// - LadderFact calculates the factorial contribution of the ratio from total boson number.
		double LadderFact(int n_new, int n_old);
		// - PsiUpdate will update the variational parameters of the wavefunction when given a vector of changes.
		void PsiUpdate(VectorXd dP);
		// - PrepPsi will load any local configuration information used in the wavefunction constituents.
		void PrepPsi(Config* Cfg);
		// - PsiCfgUpdate will update any local configuration information after a configuration change.
		void PsiCfgUpdate(); // Proposed local information changes contained within each Modifier.
		// If required later, will also have PsiCfgUpdateExt which uses an external local information source.
		// - PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
		// -- appropriate local information for one configuration loaded into the Ansatz.
		double PsiRatio(Diff diff);
		// If required later, will also have PsiRatioExt which loads local information to an external object.
		// - LogDeriv will return the logarithmic derivatives of the wavefunction w.r.t. its parameters.
		VectorXd LogDeriv(Config* Cfg) const;
	};
}

#endif