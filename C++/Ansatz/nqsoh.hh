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
		NQSOH(Hilbert* hlb, Graph* grp, int HUDen, vector<double> StartParams, vector<double> noise);
		// - Initialisation with pre-existing parameters.
		NQSOH(Hilbert* hlb, Graph* grp, int HUDen, VectorXd ParameterList);
		// Internal parameter organisation functions.
		// - NetworkInit will initialise the visible and hidden units of the RBM.
		void NetworkInit(Hilbert* hlb, Graph* grp, int HUDen);
		// - ParamFill will populate the vectors and matrices given relevant parameters.
		void ParamFill();
		// - ParamCheck will test the values of all parameters and ensure that they are less than the cap.
		void ParamCheck();
		// Observer functions:
		// - Nvisible returns the number of visible sites.
		int Nvisible() const;
		// - Nhidden returns the number of hidden sites.
		int Nhidden() const;
		// - Nparam returns the number of parameters.
		int Nparam() const;
		// - VarFlag will return 0 or 1 according to whether variational modification is permitted.
		bool VarFlag() const;
		// - OptIndRead returns the individual parameter flag vector.
		VectorXd OptIndRead() const;
		// - ParamLoad will change the existing parameters to the provided parameters.
		void ParamLoad(VectorXd NewParams);
		// - ParamList returns the parameters in the Modifier
		// -- List order is <a, b, Wv>.
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
		
		// - CfgRead will allow for additional modification of FullCfg methods in Config.
		VectorXd CfgRead(Config* Cfg) const; // Reads in configuration and outputs one-hot version.
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