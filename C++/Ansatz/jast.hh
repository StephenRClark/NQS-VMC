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
		Jast(Hilbert* hlb, Graph* grp, double J0, vector<double> noise);
		// - Initialisation with pre-existing parameters.
		Jast(Hilbert* hlb, Graph* grp, VectorXd ParameterList);
		// Internal parameter organisation functions.
		// - MatrixInit will initialise the Jastrow parameter index matrix.
		void MatrixInit(Hilbert* hlb, Graph* grp);
		// - ParamFill will populate the vectors and matrices given relevant parameters.
		void ParamFill();
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
		void ParamLoad(VectorXd NewParams);
		// - ParamList returns the parameters in the Modifier.
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