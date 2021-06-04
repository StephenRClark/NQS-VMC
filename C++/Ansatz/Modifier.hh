// ------------------------------------------------------------------------------------------------------------
// Modifier.hh - contains the definitions of the Modifier object class, which is used as a component for the
//				 Ansatz object class. While all methods are also present in Reference, Modifier is written as 
//				 a separate class to avoid stacking of multiple References and allow multiple Modifiers.
// ------------------------------------------------------------------------------------------------------------

#ifndef MODIFIER_HH
#define MODIFIER_HH

#include <Eigen/Dense>
#include "Hilbert/Hilbert.hh"

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	class Modifier // Usage of Modifier subclasses will be generating new pointers, so no implementer subclass.
	{
	public:
		// Observer functions for extracting information:
		// - Nparam will return the number of variational parameters.
		virtual int Nparam() const = 0;
		// - VarFlag will return 0 or 1 according to whether variational modification is permitted.
		virtual bool VarFlag() const = 0;
		// - OptIndRead will return individual parameter flags for variational modification.
		virtual VectorXd OptIndRead() const = 0;
		// - ParamList will return a vector of the variational parameters in the Reference.
		virtual VectorXd ParamList() const = 0;

		// Variational modification management functions:
		// - VarSwitch will alter VFlag if allowed.
		virtual void VarSwitch() = 0;
		// - OptIndLoad will load a vector of updated parameter flags (values 0/1).
		virtual void OptIndLoad(VectorXd newinds) = 0;
		// - RndBtcSelect will randomly alter the individual parameter flags.
		// -- Default zero argument version will select 1/e of the parameters.
		virtual void RndBtcSelect() = 0;
		// -- Second version with fraction will select input fraction if valid.
		virtual void RndBtcSelect(double pfrac) = 0;
		// - SetParamCap will change the maximum allowed amplitude of a single parameter.
		virtual void SetParamCap(double newcap) = 0;
		// - ParamLoad will change the existing parameters to the provided parameters.
		virtual void ParamLoad(VectorXd NewParams) = 0;

		// Wavefunction calculation functions:
		// - CfgRead will allow for additional modification of FullCfg methods in Config.
		virtual VectorXd CfgRead(Config* Cfg) const = 0;
		// - PsiUpdate will update the variational parameters of the wavefunction when given a vector of changes.
		virtual void PsiUpdate(VectorXd dP) = 0;
		// - PrepPsi will load any local configuration information used in the wavefunction constituents.
		virtual void PrepPsi(Config* Cfg) = 0;
		// - PsiCfgUpdate will update any local configuration information after a configuration change.
		virtual void PsiCfgUpdate() = 0; // Proposed local information changes contained within each Modifier.
		// If required later, will also have PsiCfgUpdateExt which uses an external local information source.
		// - PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
		// -- appropriate local information for one configuration loaded into the Ansatz.
		virtual double PsiRatio(Diff diff) = 0;
		// If required later, will also have PsiRatioExt which loads local information to an external object.
		// - LogDeriv will return the logarithmic derivatives of the wavefunction w.r.t. its parameters.
		virtual VectorXd LogDeriv(Config* Cfg) const = 0;
	};
}

#include "Ansatz/nqshdv.hh"
#include "Ansatz/nqsmh.hh"
#include "Ansatz/nqsoh.hh"
#include "Ansatz/gutz.hh"
#include "Ansatz/jast.hh"
#include "Ansatz/mbc.hh"

#endif