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
		NNMB(Hilbert* hlb, Graph* grp, double G0, vector<double> noise);
		// - Initialisation with pre-existing parameters.
		NNMB(Hilbert* hlb, Graph* grp, VectorXd ParameterList);
		// Internal parameter organisation functions.
		// - ListInit will initialise the nearest neighbour index lists.
		void ListInit(Hilbert* hlb, Graph* grp);
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