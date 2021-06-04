// ------------------------------------------------------------------------------------------------------------
// StochReconfig.hh - contains the definitions of the StochReconfig optimiser object, which utilises stochastic
//					  reconfiguration to optimise a provided Ansatz.
// ------------------------------------------------------------------------------------------------------------

#ifndef STOCHASTIC_RECONFIG_HH
#define STOCHASTIC_RECONFIG_HH

#include <iostream>
#include <vector>
#include <random>
#include <Eigen/Dense>
#include "Ansatz/Ansatz.hh"
#include "Operator/Operator.hh"
#include "Sampler/MCSampler.hh"

namespace nqsvmc
{
	// Forward declaration of the StochReconfig class.
	class StochReconfig
	{
	private:
		// General optimiser parameters:
		int Npass; // Number of optimisation passes.
		int Ncore; // Number of cores available for multithreading.
		int ExtraSamp; // Number of extra samples to take if instability is detected.
		double LearnRate; // Parameter change magnitude multiplier.
		double BFrac = exp(-1.0); // Parameter fraction for random batch optimisation.
		// Tolerances in energy and parameter changes:
		double EMTol = 5; // Maximum energy magnitude tolerance. Adjust as the physics requires.
		double dEVTol = 0.25; // Maximum tolerance in absolute value energy increase.
		double dERTol = 0.05; // Maximum tolerance in relative energy increase.
		double dEMTol = 1e-6; // Tolerance in energy change to determine convergence.
		double dPTol = 1e-10; // Minimum parameter change tolerance.
		double PVTol = 1e-10; // Minimum parameter value tolerance.
		double MRTol = 0.0001; // Tolerance in move acceptance ratio.
		// Internal flags counters for error checking:
		int itcount = 0; // Counts iterations in extra sampling phase.
		int encount = 0; // Counts iterations with energy variation below threshold.
		int mrcount = 0; // Counts iterations with move acceptance rate below threshold.
		bool SampFlag = 0; // Activates extra sampling phase.
		bool TermFlag = 0; // Terminates runs early under certain conditions.
		// Stochastic reconfiguration regularisation parameters.
		double Rmax; // Maximum or starting regularisation strength.
		double Rmin; // Minimum regularisation strength.
		double Rmult; // Regularisation attenuation multiplier.
		double STol = 1e-20; // Tolerance in overlap matrix element magnitudes.
		int PSave = 10; // Number of steps backward to take in the event of an error.
		int PShift = 10; // Number of multiplier factors to roll back in the event of an error.
	public:
		// Constructors for the StochReconfig object:
		// - Basic variant with only number of passes and number of cores to utilise.
		StochReconfig(int NumPass, int NumCore);
		// - More advanced variant that specifies extra samples, learn rate and regularisation.
		StochReconfig(int NumPass, int NumCore, int ExtraSamples, double LearningRate, vector<double> RegVec);
		// - Complete definition initialisation - used when loading in from external.
		StochReconfig(int NumPass, int NumCore, int ExtraSamples, double LearningRate, vector<double> RegVec,
			vector<int> RollbackVals, vector<double> EneValTolerances, vector<double> EneChgTolerances, vector<double> ParamTolerances);
		// Property management functions:
		// - SetNpass assigns the number of optimisation passes to apply to the Ansatz.
		void SetNpass(int Nnew);
		// - SetExtraSamples assigns the number of additional samples to use after a failed pass.
		void SetExtraSamples(int Nextra);
		// - SetEVTolerances assigns the maximum magnitude and minimum change in energy tolerances.
		void SetEVTolerances(double EMagTol, double dEMagTol);
		// - SetdETolerances assigns the relative and absolute positive energy change tolerances.
		void SetdETolerances(double RelTol, double AbsTol);
		// - SetParamTolerances assigns the tolerance in parameter change and total magnitudes.
		void SetParamTolerances(double MagTol, double ChgTol);
		// - SetBatchFraction assigns the fraction of parameters for random batch selection.
		void SetBatchFraction(double BatchFrac);
		// - SetRegularisation assigns the regularisation terms Rmax, Rmin and Rmult.
		void SetRegularisation(double R_maximum, double R_minimum, double R_multiplier);
		// - SetLearnRate will alter the parameter step size factor.
		void SetLearnRate(double L_new);
		// - SetRollback will alter the rollback parameters PSave and PShift.
		void SetRollback(int psave, int pshift);

		// Calculation functions for matrix checking, ansatz rollback and parameter calculation.
		// - SMatCheck will mark out elements of the S matrix below a magnitude threshold.
		VectorXi SMatCheck(MatrixXd S_in) const;
		// - ParamCalc takes the S matrix, F vector and a list of indices to output parameter changes dP.
		VectorXd ParamCalc(MatrixXd Smat, VectorXd Fvec, VectorXi S_inds) const;
		// - AnsRollback will select the parameter list that has given the lowest energy and load them into the Ansatz.
		void AnsRollback(Ansatz* AnsObj, vector<double> EneLog, vector<VectorXd> ParamLog) const;
		// - EnergyCheck will try to identify any outrageous changes in energy and correct them.
		void EnergyCheck(double E_p, double E_ref, Ansatz* AnsObj, int p,
			double MRate, int &Pshift, vector<double> EneLog, vector<VectorXd> ParamLog);
		// - SampleReset will simply zero out variables used to store sampled quantities.
		void SampleReset(double& EnAvg, VectorXd& dLogp, double& MRate, vector<MatrixXd>& OpVals) const;
		// - FlagReset sets all internal counters to zero / false.
		void FlagReset();

		// Optimisation functions:
		// - Optimise will alter the provided Ansatz using parameter changes calculated.
		VectorXd Optimise(Ansatz* AnsObj, MCSampler* SamplerObj);
		// - RndBtcOptimise will select a random batch of parameters to optimise at each step.
		VectorXd RndBtcOptimise(Ansatz* AnsObj, MCSampler* SamplerObj);
		// - PrmAvgOptimise will average out the parameter values calculated over all the passes.
		VectorXd PrmAvgOptimise(Ansatz* AnsObj, MCSampler* SamplerObj);		
	};
}

#endif