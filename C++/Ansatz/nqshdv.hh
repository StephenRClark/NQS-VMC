// ------------------------------------------------------------------------------------------------------------
// nqshdv.hh - contains the definitions of the NQS class with multinomial hidden units, which is a subclass of 
//			   Modifier. This variant will assume translation invariance, and is intended for use with bosonic
//			   systems.
// ------------------------------------------------------------------------------------------------------------

#ifndef NQS_HDV_HH
#define NQS_HDV_HH

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
	class NQSSH : public Modifier
	{
	private: // Properties associated with the spin-like hidden unit NQS.
		Graph* g_; // Pointer to a Graph object for translation invariance.
		Hilbert* h_; // Pointer to a Hilbert object.
		// Network parameters:
		int Nv; // Number of visible sites, scalar.
		int Nh; // Number of hidden units, scalar.
		int Alpha; // Hidden unit density / unique stack number, scalar.
		int HDim; // Hidden unit dimension, scalar.
		// Variational parameters:
		double a; // Linear visible bias parameter, scalar.
		double A; // Non-linear visible bias parameter, scalar.
		VectorXd b; // Linear hidden bias parameters, Alpha x 1 vector.
		VectorXd B; // Non-linear hidden bias parameters, Alpha x 1 vector.
		MatrixXd Wv; // Hidden-visible coupling parameters, Alpha x Nv matrix.
		// Full networked parameter lists:
		VectorXd av; // Linear visible bias, Nv x 1 vector.
		VectorXd Av; // Non-linear visible bias, Nv x 1 vector.		
		VectorXd bv; // Linear hidden bias, Nh x 1 vector.
		VectorXd Bv; // Non-linear hidden bias, Nh x 1 vector.		
		MatrixXd Wm; // Hidden-visible couplings, Nh x Nv matrix.
		// Local information:
		VectorXd Nsq; // Vector of squared visible occupancies, Nv x 1 vector.
		VectorXd Theta; // Effective angles, Nh x 1 vector.
		VectorXd NsqP; // New proposed squared visible occupancies.
		VectorXd ThetaP; // New proposed effective angles.
		// Variational management:
		double ParamCap; // Maximum parameter magnitude, scalar.
		int Np; // Number of parameters.
		bool VFlag; // Variational modification activation flag.
		VectorXd OptInds; // Individual parameter flags, Np x 1 Boolean vector.
	public:
		// Constructors for the spin-hidden NQS:
		// - Random initialisation with starting values of each parameter type.
		NQSSH(Hilbert* hlb, Graph* grp, int HUDen, vector<double> StartParams, vector<double> noise);
		// - Initialisation with pre-existing parameters.
		NQSSH(Hilbert* hlb, Graph* grp, int HUDen, VectorXd ParameterList);
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
		// Generalised hidden unit trace functions:
		VectorXd SHTrace(VectorXd theta_in) const; // Calculate trace over hidden units for given Theta.
		VectorXd dT_SHTrace(VectorXd theta_in) const; // Calculate derivative of trace w.r.t Theta.
		VectorXd dB_SHTrace(VectorXd theta_in) const; // Calculate derivative of trace w.r.t B.
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

	class NQSNH : public Modifier
	{
	private: // Properties associated with the spin-like hidden unit NQS.
		Graph* g_; // Pointer to a Graph object for translation invariance.
		Hilbert* h_; // Pointer to a Hilbert object.
		// Network parameters:
		int Nv; // Number of visible sites, scalar.
		int Nh; // Number of hidden units, scalar.
		int Alpha; // Hidden unit density / unique stack number, scalar.
		int HDim; // Hidden unit dimension, scalar.
		// Variational parameters:
		double a; // Linear visible bias parameter, scalar.
		double A; // Non-linear visible bias parameter, scalar.
		VectorXd b; // Linear hidden bias parameters, Alpha x 1 vector.
		VectorXd B; // Non-linear hidden bias parameters, Alpha x 1 vector.
		MatrixXd Wv; // Hidden-visible coupling parameters, Alpha x Nv matrix.
		// Full networked parameter lists:
		VectorXd av; // Linear visible bias, Nv x 1 vector.
		VectorXd Av; // Non-linear visible bias, Nv x 1 vector.		
		VectorXd bv; // Linear hidden bias, Nh x 1 vector.
		VectorXd Bv; // Non-linear hidden bias, Nh x 1 vector.		
		MatrixXd Wm; // Hidden-visible couplings, Nh x Nv matrix.
		// Local information:
		VectorXd Nsq; // Vector of squared visible occupancies, Nv x 1 vector.
		VectorXd Theta; // Effective angles, Nh x 1 vector.
		VectorXd NsqP; // New proposed squared visible occupancies.
		VectorXd ThetaP; // New proposed effective angles.
		// Variational management:
		double ParamCap; // Maximum parameter magnitude, scalar.
		int Np; // Number of parameters.
		bool VFlag; // Variational modification activation flag.
		VectorXd OptInds; // Individual parameter flags, Np x 1 Boolean vector.
	public:
		// Constructors for the number-hidden NQS:
		// - Random initialisation with starting values of each parameter type.
		NQSNH(Hilbert* hlb, Graph* grp, int HUDen, vector<double> StartParams, vector<double> noise);
		// - Initialisation with pre-existing parameters.
		NQSNH(Hilbert* hlb, Graph* grp, int HUDen, VectorXd ParameterList);
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
		// Generalised hidden unit trace functions:
		VectorXd NHTrace(VectorXd theta_in) const; // Calculate trace over hidden units for given Theta.
		VectorXd dT_NHTrace(VectorXd theta_in) const; // Calculate derivative of trace w.r.t Theta.
		VectorXd dB_NHTrace(VectorXd theta_in) const; // Calculate derivative of trace w.r.t B.
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