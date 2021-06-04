// ------------------------------------------------------------------------------------------------------------
// Ansatz.hh - contains the declarations of the Ansatz object class, outlining its interactions with the Config
//			   type of object. The Ansatz class contains all the methods of calculating wavefunction 
//			   properties such as amplitudes and parameter derivatives. Within the Ansatz is a Reference 
//			   object and some number of Modifier objects, and the methods in Ansatz will call the methods of 
//			   its constituents.
// ------------------------------------------------------------------------------------------------------------

#ifndef ANSATZ_HH
#define ANSATZ_HH

#include <iostream>
#include <Eigen/Dense>
#include "Hilbert/Hilbert.hh"
#include "Reference.hh"
#include "Modifier.hh"

namespace nqsvmc
{
	// Full declaration of functions in Ansatz to allow other parts of the code to interface with Ansatz.
	class Ansatz
	{
	private:
		Reference* r_; // Pointer to Reference object.
		vector<Modifier*> m_; // Vector of pointers to Modifier objects.
		Hilbert* h_; // Pointer to Hilbert object used to check compatibility.
		int Nvar; // Number of variational components (Reference/Modifier).
		int NpTotal; // Total number of variational parameters.
	public:
		// Constructors for the Ansatz object class
		// - Constructor with inputs of pointers.
		Ansatz(Reference* Ref, vector<Modifier*> Mods, Hilbert* Hlb);

		// Observer functions for extracting information from an Ansatz object:
		// - NumVar returns the number of variational components in the Ansatz.
		int NumVar() const;
		// - Nparam returns the number of variational parameters in the Ansatz.
		int Nparam() const;
		// - Nmod returns the number of Modifiers.
		int Nmod() const;
		// - ParamList will return a vector of parameters in the ansatz from References and Modifiers.
		VectorXd ParamList() const;

		// Variational management functions:
		// - RndBtcSelect will randomly select parameters within the Ansatz components.
		// -- First version with no specified batch fraction.
		void RndBtcSelect();
		// -- Second version uses specified batch fraction.
		void RndBtcSelect(double pfrac);

		// Wavefunction calculation and update functions:
		// - PsiUpdate will update the variational parameters of the wavefunction when given a vector of changes.
		void PsiUpdate(VectorXd dP);
		// - ParamLoad will overwrite the variatonal parameters of the wavefunction when given a vector of parameters.
		void ParamLoad(VectorXd P);
		// - OptIndLoad will overwrite the variatonal parameter flags of the wavefunction when given a vector of new indices.
		void OptIndLoad(VectorXd newinds);
		// - PrepPsi will load any local configuration information used in the wavefunction constituents.
		void PrepPsi(Config* Cfg);
		// - PsiCfgUpdate will update any local configuration information after a configuration change.
		void PsiCfgUpdate();
		// Proposed local information changes contained within each Modifier.
		// If required later, will also have PsiCfgUpdateExt which uses an external local information source.
		// - PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
		// -- appropriate local information for one configuration loaded into the Ansatz.
		double PsiRatio(Diff diff);
		// If required later, will also have PsiRatioExt which loads local information to an external object.
		// - LogDeriv will return the logarithmic derivatives of the wavefunction w.r.t. its parameters.
		VectorXd LogDeriv(Config* Cfg);
	};	
}

#endif