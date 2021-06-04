// ------------------------------------------------------------------------------------------------------------
// StochReconfig.cpp - contains the definitions of the StochReconfig optimiser object, which utilises 
//					   stochastic reconfiguration to optimise a provided Ansatz.
// ------------------------------------------------------------------------------------------------------------

#ifndef STOCHASTIC_RECONFIG_CPP
#define STOCHASTIC_RECONFIG_CPP

#include <iostream>
#include <vector>
#include <random>
#include <Eigen/Dense>
#include "StochReconfig.hh"

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	// - Basic variant with only number of passes and number of cores to utilise.
	StochReconfig::StochReconfig(int NumPass, int NumCore)
	{
		SetNpass(NumPass);
		if (NumCore <= 0)
		{
			cerr << "Number of cores must be a positive integer." << endl;
			std::abort();
		}
		Ncore = NumCore;
		ExtraSamp = 5000;
		LearnRate = 0.1;
		Rmax = 10000;
		Rmin = 0.0001;
		Rmult = 0.9;
		return;
	}
	// - More advanced variant that specifies extra samples, learn rate and regularisation.
	StochReconfig::StochReconfig(int NumPass, int NumCore, int ExtraSamples, double LearningRate, vector<double> RegVec)
	{
		SetNpass(NumPass);
		if (NumCore <= 0)
		{
			cerr << "Number of cores must be a positive integer." << endl;
			std::abort();
		}
		Ncore = NumCore;
		SetExtraSamples(ExtraSamples);
		SetLearnRate(LearningRate);
		if (RegVec.size() != 3)
		{
			cerr << "Vector of regularisation parameters requires three terms - maximum, minumum and multiplier." << endl;
			std::abort();
		}
		SetRegularisation(RegVec[0], RegVec[1], RegVec[2]);
		return;
	}
	// - Complete definition initialisation - used when loading in from external.
	StochReconfig::StochReconfig(int NumPass, int NumCore, int ExtraSamples, double LearningRate, vector<double> RegVec,
		vector<int> RollbackVals, vector<double> EneValTolerances, vector<double> EneChgTolerances, vector<double> ParamTolerances)
	{
		SetNpass(NumPass);
		if (NumCore <= 0)
		{
			cerr << "Number of cores must be a positive integer." << endl;
			std::abort();
		}
		Ncore = NumCore;
		SetExtraSamples(ExtraSamples);
		SetLearnRate(LearningRate);
		if (RegVec.size() != 3)
		{
			cerr << "Vector of regularisation parameters requires three terms - maximum, minumum and multiplier." << endl;
			std::abort();
		}
		SetRegularisation(RegVec[0], RegVec[1], RegVec[2]);
		if (RollbackVals.size() != 2)
		{
			cerr << "Vector of rollback parameters requires two terms - number of saved parameter sets and regularisation rollback factor." << endl;
			std::abort();
		}
		SetRollback(RollbackVals[0], RollbackVals[1]);		
		if (EneValTolerances.size() != 2)
		{
			cerr << "Vector of energy value tolerances requires two terms - absolute value tolerance and minimum change tolerance." << endl;
			std::abort();
		}
		SetEVTolerances(EneValTolerances[0], EneValTolerances[1]);
		if (EneChgTolerances.size() != 2)
		{
			cerr << "Vector of energy change tolerances requires two terms - maximum relative change and maximum value change." << endl;
			std::abort();
		}
		SetdETolerances(EneChgTolerances[0], EneChgTolerances[1]);
		if (ParamTolerances.size() != 2)
		{
			cerr << "Vector of parameter tolerances requires two terms - minimum parameter magnitude and minimum parameter change magnitude." << endl;
			std::abort();
		}
		SetParamTolerances(ParamTolerances[0], ParamTolerances[1]);
		return;
	}

	// Property management functions:
	// - SetNpass assigns the number of optimisation passes to apply to the Ansatz.
	void StochReconfig::SetNpass(int Nnew)
	{
		if (Nnew <= 0)
		{
			cerr << "New number of passes must be a positive integer." << endl;
			std::abort();
		}
		Npass = Nnew;
		return;
	}
	// - SetExtraSamples assigns the number of additional samples to use after a failed pass.
	void StochReconfig::SetExtraSamples(int Nextra)
	{
		if (Nextra < 0)
		{
			cerr << "New number of extra samples must be a non-negative integer." << endl;
			std::abort();
		}
		ExtraSamp = Nextra;
		return;
	}
	// - SetEVTolerances assigns the maximum magnitude and minimum change in energy tolerances.
	void StochReconfig::SetEVTolerances(double EMagTol, double dEMagTol)
	{
		if (EMagTol < 0)
		{
			cout << "Taking absolute value of the energy magnitude tolerance input..." << endl;
			EMagTol = abs(EMagTol);
		}
		if (dEMagTol < 0)
		{
			cout << "Taking absolute value of the minimum energy change tolerance input..." << endl;
			dEMagTol = abs(dEMagTol);
		}
		EMTol = EMagTol;
		dEMTol = dEMagTol;
		return;
	}
	// - SetdETolerances assigns the relative and absolute positive energy change tolerances.
	void StochReconfig::SetdETolerances(double RelTol, double AbsTol)
	{
		if (RelTol < 0)
		{
			cout << "Taking absolute value of the relative tolerance input..." << endl;
			RelTol = abs(RelTol);
		}
		if (AbsTol < 0)
		{
			cout << "Taking absolute value of the absolute tolerance input..." << endl;
			AbsTol = abs(AbsTol);
		}
		dERTol = RelTol;
		dEVTol = AbsTol;
		return;
	}
	// - SetParamTolerances assigns the tolerance in parameter change and total magnitudes.
	void StochReconfig::SetParamTolerances(double MagTol, double ChgTol)
	{
		if (MagTol < 0)
		{
			cout << "Taking absolute value of the parameter magnitude tolerance input..." << endl;
			MagTol = abs(MagTol);
		}
		if (ChgTol < 0)
		{
			cout << "Taking absolute value of the parameter change magnitude tolerance input..." << endl;
			ChgTol = abs(ChgTol);
		}
		dPTol = ChgTol;
		PVTol = MagTol;
		return;
	}
	// - SetBatchFraction assigns the fraction of parameters for random batch selection.
	void StochReconfig::SetBatchFraction(double BatchFrac)
	{
		if ((BatchFrac <= 0) || BatchFrac > 1)
		{
			cerr << "Batch fraction must be a positive number between 0 and 1 (1 inclusive)." << endl;
			std::abort();
		}
		BFrac = BatchFrac;
		return;
	}
	// - SetRegularisation assigns the regularisation terms Rmax, Rmin and Rmult.
	void StochReconfig::SetRegularisation(double R_maximum, double R_minimum, double R_multiplier)
	{
		if (R_maximum < 0)
		{
			cout << "Taking absolute value of the maximum regularisation input..." << endl;
			R_maximum = abs(R_maximum);
		}
		if (R_minimum < 0)
		{
			cout << "Taking absolute value of the minimum regularisation input..." << endl;
			R_minimum = abs(R_minimum);
		}
		if (R_multiplier < 0)
		{
			cout << "Taking absolute value of the regularisation multiplier input..." << endl;
			R_multiplier = abs(R_multiplier);
		}
		if (R_multiplier == 0)
		{
			cerr << "Multiplier cannot be zero." << endl;
			std::abort();
		}
		Rmax = R_maximum;
		Rmin = R_minimum;
		Rmult = R_multiplier;
		return;
	}
	// - SetLearnRate will alter the parameter step size factor.
	void StochReconfig::SetLearnRate(double L_new)
	{
		if (L_new <= 0)
		{
			cerr << "New learn rate must be a real positive number." << endl;
			std::abort();
		}
		LearnRate = L_new;
		return;
	}
	// - SetRollback will alter the rollback parameters PSave and PShift.
	void StochReconfig::SetRollback(int psave, int pshift)
	{
		if (psave < 5)
		{
			std::cout << "PSave minimum is 5 passes. Setting to 5." << endl;
			PSave = 5;
		}
		else
		{
			PSave = psave;
		}
		PShift = pshift;
		return;
	}

	// Calculation functions for matrix checking, ansatz rollback and parameter calculation.
	// - SMatCheck will mark out elements of the S matrix below a magnitude threshold.
	VectorXi StochReconfig::SMatCheck(MatrixXd S_in) const
	{
		// Judgement will be made on the values of the diagonal of S_in.
		VectorXd S_diag = S_in.diagonal();
		VectorXi S_inds = VectorXi::Zero(S_diag.size());
		for (int s = 0; s < S_diag.size(); s++)
		{
			if (S_diag(s) > STol)
			{
				S_inds(s) = 1;
			}
		}
		return S_inds;
	}

	// - ParamCalc takes the S matrix, F vector and a list of indices to output parameter changes dP.
	VectorXd StochReconfig::ParamCalc(MatrixXd Smat, VectorXd Fvec, VectorXi S_inds) const
	{
		int ns = S_inds.sum();
		int Nparam = (int)S_inds.size();
		MatrixXd S_fin = MatrixXd::Zero(ns, ns);
		VectorXd F_fin = VectorXd::Zero(ns);
		int nc = 0;
		for (int n = 0; n < Nparam; n++)
		{
			int mc = 0;
			if (S_inds(n) > 0)
			{
				F_fin(nc) = Fvec(n);
				for (int m = 0; m < Nparam; m++)
				{
					if (S_inds(m) > 0)
					{
						S_fin(nc, mc) = Smat(n, m);
						mc++;
					}
				}
				nc++;
			}
		}
		S_fin = 0.5 * (S_fin + S_fin.transpose());
		// VectorXd S_diag = S_fin.diagonal().array().sqrt();
		// MatrixXd S_scale = S_diag * S_diag.transpose();
		// F_fin = F_fin.array() / S_diag.array();
		// S_fin = S_fin.array() / S_scale.array();
		VectorXd dP_raw = (S_fin.bdcSvd(ComputeThinU | ComputeThinV).solve(F_fin));
		// VectorXd dP_raw = (S_fin.fullPivLu().solve(F_fin));
		// dP_raw = dP_raw.array() / S_diag.array();
		VectorXd dP_fin = VectorXd::Zero(Nparam);
		int pc = 0;
		for (int p = 0; p < Nparam; p++)
		{
			if (S_inds(p) > 0)
			{
				if ((isnan(dP_raw(pc)) == 0) || (isinf(dP_raw(pc)) == 0))
				{
					// Only add to final vector if not NaN or Inf.
					dP_fin(p) = dP_raw(pc);
				}
				pc++;
			}
		}
		dP_fin = -dP_fin.array() * LearnRate;
		return dP_fin;
	}
	// - AnsRollback will select the parameter list that has given the lowest energy and load them into the Ansatz.
	void StochReconfig::AnsRollback(Ansatz* AnsObj, vector<double> EneLog, vector<VectorXd> ParamLog) const
	{
		int P_ind = (int)std::distance(EneLog.begin(), std::min_element(EneLog.begin(), EneLog.end()));
		cout << "Rolling back to entry " << P_ind << " in logs, energy " << EneLog[P_ind] << ".";
		VectorXd Params = ParamLog[P_ind];
		VectorXd OptInds = VectorXd::Zero(Params.size());
		for (int p = 0; p < Params.size(); p++)
		{
			if (abs(Params(p)) > 1e-30)
			{
				Params(p) = 0;
				OptInds(p) = 1;
			}
		}
		AnsObj->ParamLoad(ParamLog[P_ind]);
		AnsObj->OptIndLoad(OptInds);
		return;
	}
	// - EnergyCheck will try to identify any outrageous changes in energy and correct them.
	void StochReconfig::EnergyCheck(double E_p, double E_ref, Ansatz* AnsObj, int p,
		double MRate, int &Pshift, vector<double> EneLog, vector<VectorXd> ParamLog)
	{
		double dEV = E_p - E_ref;
		double dER = dEV / abs(EneLog.back());
		if ((abs(dER) < dEMTol) && (p > (Npass / 5)))
		{
			encount++;
		}
		else
		{
			encount = 0;
		}
		if (SampFlag)
		{
			itcount++;
			if (itcount == PSave)
			{
				itcount = 0;
				cout << "Ending extra sampling period." << endl;
				SampFlag = false;
			}
		}
		if (MRate < MRTol)
		{
			mrcount++;
		}
		else
		{
			mrcount = 0;
		}
		if (SampFlag)
		{
			if ((dEV > dEVTol) && (dER > dERTol))
			{
				cout << "Unfavourable energy proposed during extra sampling - backtracking." << endl;
				itcount = 0;
				AnsRollback(AnsObj, EneLog, ParamLog);
			}
			else if (mrcount > PSave)
			{
				cout << "Move acceptance rate has fallen below threshold for " << PSave <<
					" passes during extra sampling - backtracking." << endl;
				itcount = 0;
				AnsRollback(AnsObj, EneLog, ParamLog);
			}
			else if ((abs(dER) < dEMTol) && (encount == PSave))
			{
				cout << "Energy change has fallen below threshold for " << PSave <<
					"consecutive runs - assuming stuck in minimum and backtracking." << endl;
				itcount = 0;
				AnsRollback(AnsObj, EneLog, ParamLog);
			}
		}
		else
		{
			if ((dEV > dEVTol) && (dER > dERTol))
			{
				cout << "Unfavourable energy proposed - increasing sampling and backtracking." << endl;
				SampFlag = true; Pshift += PShift;
				itcount = 0;
				AnsRollback(AnsObj, EneLog, ParamLog);
			}
			else if (mrcount > PSave)
			{
				cout << "Move acceptance rate has fallen below threshold for " << PSave <<
					" passes - increasing sampling and backtracking." << endl;
				SampFlag = true; Pshift += PShift;
				itcount = 0;
				AnsRollback(AnsObj, EneLog, ParamLog);
			}
			else if ((abs(dER) < dEMTol) && (encount == PSave))
			{
				cout << "Energy change has fallen below threshold for " << PSave <<
					" consecutive passes - assuming convergence has been reached and terminating further passes." << endl;
				TermFlag = true;
			}
		}
		return;
	}
	void StochReconfig::SampleReset(double& EnAvg, VectorXd& dLogp, double& MRate, vector<MatrixXd>& OpVals) const
	{
		EnAvg *= 0;
		dLogp.setZero();
		MRate *= 0;
		for (int o = 0; o < OpVals.size(); o++)
		{
			OpVals[o].setZero();
		}
		return;
	}
	// - FlagReset sets all internal counters to zero / false.
	void StochReconfig::FlagReset()
	{
		SampFlag = 0;
		TermFlag = 0;
		itcount = 0;
		encount = 0;
		mrcount = 0;
		return;
	}

	// Optimisation functions:
	// - Optimise will alter the provided Ansatz using parameter changes calculated.
	VectorXd StochReconfig::Optimise(Ansatz* AnsObj, MCSampler* SamplerObj)
	{
		int Nparam = AnsObj->Nparam();
		Operator EdLogp = Operator("EdP");
		Operator dLogpq = Operator("dPdQ");
		vector<Operator*> grad_op{ &EdLogp,&dLogpq };
		MatrixXd EdP_val = MatrixXd::Zero(Nparam, 1);
		MatrixXd dPdQ_val = MatrixXd::Zero(Nparam, Nparam);
		vector<MatrixXd> grad_val{ EdP_val,dPdQ_val };
		// Preparing storage and error counters.
		VectorXd EnIter = VectorXd::Zero(Npass);
		if (PSave < 5)
		{
			PSave = 5;
		}
		vector<VectorXd> ParamLog(PSave, EdP_val);
		vector<double> EneLog(PSave);
		// Reset internal flags and counters.
		FlagReset();
		int Nsamp0 = SamplerObj->NumSamp();
		int NsampP = Nsamp0;
		int Pshift = 0;
		double EnP = 0;
		double MRate = 0;
		double Rp = 0;
		VectorXd dLogp = VectorXd::Zero(Nparam);
		VectorXd dP = VectorXd::Zero(Nparam);
		for (int p = 0; p < Npass; p++)
		{
			NsampP = Nsamp0 + SampFlag * ExtraSamp;
			// Reset values each loop.
			SampleReset(EnP, dLogp, MRate, grad_val);
			SamplerObj->MCMCSample(AnsObj, EnP, dLogp, grad_op, grad_val, MRate);
			// Construct S matrix and F vector.
			MatrixXd Smat_raw = grad_val[1] - (dLogp * dLogp.transpose());
			// Apply matrix check to trim zero / small entries later.
			VectorXi S_inds = SMatCheck(Smat_raw);
			VectorXd Smat_diag = Smat_raw.diagonal();
			VectorXd Fvec_raw = grad_val[0] - (EnP * dLogp);			
			// Calculate regularisation.
			Rp = max(Rmin, (Rmax * pow(Rmult, p + Pshift)));
			MatrixXd Sdiag_reg = Rp * Smat_diag.asDiagonal();
			MatrixXd Smat_reg = Smat_raw + Sdiag_reg;			
			// Use S_inds to reconstruct the S matrix without omitted entries.
			int ns = S_inds.sum();
			if (ns > 0)
			{
				dP = ParamCalc(Smat_reg, Fvec_raw, S_inds);
			}
			else
			{
				cout << "All S matrix entries below threshold - no parameter changes applied." << endl;
				dP = VectorXd::Zero(Nparam);
			}
			if (isnan(EnP) || isinf(EnP) || (EnP > EMTol))
			{
				cout << "Unphysical energy from sampling - proposed parameter changes rejected." << endl;
				EnP = EMTol;
				dP *= 0;
			}
			EnIter(p) = EnP;
			cout << "Energy at step " << p << ": " << EnP << endl;
			AnsObj->PsiUpdate(dP);
			if (p < PSave) // Build a log of parameters.
			{
				
				EneLog[p] = EnP;
				ParamLog[p] = AnsObj->ParamList();
			}
			else // Add to logs and perform all the relevant error checks.
			{
				EnergyCheck(EnP, EneLog[(PSave-1)], AnsObj, p, MRate, Pshift, EneLog, ParamLog);
				EneLog.erase(EneLog.begin());
				EneLog.push_back(EnP);
				ParamLog.erase(ParamLog.begin());
				ParamLog.push_back(AnsObj->ParamList());				
			}
			if (TermFlag)
			{
				EnIter = EnIter.head(p);
				break;
			}
		}
		return EnIter;
	}
	// - RndBtcOptimise will select a random batch of parameters to optimise at each step.
	VectorXd StochReconfig::RndBtcOptimise(Ansatz* AnsObj, MCSampler* SamplerObj)
	{
		int Nparam = AnsObj->Nparam();
		Operator EdLogp = Operator("EdP");
		Operator dLogpq = Operator("dPdQ");
		vector<Operator*> grad_op{ &EdLogp,&dLogpq };
		MatrixXd EdP_val = MatrixXd::Zero(Nparam, 1);
		MatrixXd dPdQ_val = MatrixXd::Zero(Nparam, Nparam);
		vector<MatrixXd> grad_val{ EdP_val,dPdQ_val };
		// Preparing storage and error counters.
		VectorXd EnIter = VectorXd::Zero(Npass);
		if (PSave < 5)
		{
			PSave = 5;
		}
		vector<VectorXd> ParamLog(PSave, EdP_val);
		vector<double> EneLog(PSave);
		// Reset internal flags and counters.
		FlagReset();
		int Nsamp0 = SamplerObj->NumSamp();
		int NsampP = Nsamp0;
		int Pshift = 0;
		double EnP = 0;
		double MRate = 0;
		double Rp = 0;
		VectorXd dLogp = VectorXd::Zero(Nparam);
		VectorXd dP = VectorXd::Zero(Nparam);
		for (int p = 0; p < Npass; p++)
		{
			NsampP = Nsamp0 + SampFlag * ExtraSamp;
			// Reset values each loop.
			SampleReset(EnP, dLogp, MRate, grad_val);
			AnsObj->RndBtcSelect(BFrac);
			SamplerObj->MCMCSample(AnsObj, EnP, dLogp, grad_op, grad_val, MRate);
			// Construct S matrix and F vector.
			MatrixXd Smat_raw = dPdQ_val - (dLogp * dLogp.transpose());
			VectorXd Smat_diag = Smat_raw.diagonal();
			VectorXd Fvec_raw = EdP_val - (EnP * dLogp);
			// Calculate regularisation.
			Rp = max(Rmin, (Rmax * pow(Rmult, p + Pshift)));
			MatrixXd Sdiag_reg = Rp * Smat_diag.asDiagonal();
			MatrixXd Smat_reg = Smat_raw + Sdiag_reg;
			// Apply matrix check to trim zero / small entries later.
			VectorXi S_inds = SMatCheck(Smat_reg);
			// Use S_inds to reconstruct the S matrix without omitted entries.
			int ns = S_inds.sum();
			if (ns > 0)
			{
				dP = ParamCalc(Smat_reg, Fvec_raw, S_inds);
			}
			else
			{
				cout << "All S matrix entries below threshold - no parameter changes applied." << endl;
				dP = VectorXd::Zero(Nparam);
			}
			if (isnan(EnP) || isinf(EnP) || (EnP > EMTol))
			{
				cout << "Unphysical energy from sampling - proposed parameter changes rejected." << endl;
				EnP = EMTol;
				dP *= 0;
			}
			EnIter(p) = EnP;
			AnsObj->PsiUpdate(dP);
			if (p < PSave) // Build a log of parameters.
			{
				EneLog[p] = EnP;
				ParamLog[p] = AnsObj->ParamList();
			}
			else // Add to logs and perform all the relevant error checks.
			{
				EnergyCheck(EnP, EneLog[(PSave-1)], AnsObj, p, MRate, Pshift, EneLog, ParamLog);
				EneLog.erase(EneLog.begin());
				EneLog.push_back(EnP);
				ParamLog.erase(ParamLog.begin());
				ParamLog.push_back(AnsObj->ParamList());				
			}
			if (TermFlag)
			{
				EnIter = EnIter.head(p);
				break;
			}
		}
		return EnIter;
	}
	// - PrmAvgOptimise will average out the parameter values calculated over all the passes.
	VectorXd StochReconfig::PrmAvgOptimise(Ansatz* AnsObj, MCSampler* SamplerObj)
	{
		int Nparam = AnsObj->Nparam();
		Operator EdLogp = Operator("EdP");
		Operator dLogpq = Operator("dPdQ");
		vector<Operator*> grad_op{ &EdLogp,&dLogpq };
		MatrixXd EdP_val = MatrixXd::Zero(Nparam, 1);
		MatrixXd dPdQ_val = MatrixXd::Zero(Nparam, Nparam);
		vector<MatrixXd> grad_val{ EdP_val,dPdQ_val };
		// Preparing storage and error counters.
		VectorXd EnIter = VectorXd::Zero(Npass);
		if (PSave < 5)
		{
			PSave = 5;
		}
		vector<VectorXd> ParamLog(PSave, EdP_val);
		vector<double> EneLog(PSave);
		// Reset internal flags and counters.
		FlagReset();
		int Nsamp0 = SamplerObj->NumSamp();
		int NsampP = Nsamp0;
		int Pshift = 0;
		double EnP = 0;
		double MRate = 0;
		double Rp = 0;
		int pcount = 0;
		VectorXd dLogp = VectorXd::Zero(Nparam);
		VectorXd P_fin = VectorXd::Zero(Nparam);
		VectorXd dP = VectorXd::Zero(Nparam);
		for (int p = 0; p < Npass; p++)
		{
			NsampP = Nsamp0 + SampFlag * ExtraSamp;
			// Reset values each loop.
			SampleReset(EnP, dLogp, MRate, grad_val);
			SamplerObj->MCMCSample(AnsObj, EnP, dLogp, grad_op, grad_val, MRate);
			// Construct S matrix and F vector.
			MatrixXd Smat_raw = dPdQ_val - (dLogp * dLogp.transpose());
			VectorXd Smat_diag = Smat_raw.diagonal();
			VectorXd Fvec_raw = EdP_val - (EnP * dLogp);
			// Calculate regularisation.
			Rp = max(Rmin, (Rmax * pow(Rmult, p + Pshift)));
			MatrixXd Sdiag_reg = Rp * Smat_diag.asDiagonal();
			MatrixXd Smat_reg = Smat_raw + Sdiag_reg;
			// Apply matrix check to trim zero / small entries later.
			VectorXi S_inds = SMatCheck(Smat_reg);
			// Use S_inds to reconstruct the S matrix without omitted entries.
			int ns = S_inds.sum();
			if (ns > 0)
			{
				dP = ParamCalc(Smat_reg, Fvec_raw, S_inds);
			}
			else
			{
				cout << "All S matrix entries below threshold - no parameter changes applied." << endl;
				dP = VectorXd::Zero(Nparam);
			}
			if (isnan(EnP) || isinf(EnP) || (EnP > EMTol))
			{
				cout << "Unphysical energy from sampling - proposed parameter changes rejected." << endl;
				EnP = EMTol;
				dP *= 0;
			}
			EnIter(p) = EnP;
			AnsObj->PsiUpdate(dP);
			if (p < PSave) // Build a log of parameters.
			{
				EneLog[p] = EnP;
				ParamLog[p] = AnsObj->ParamList();
			}
			else // Add to logs and perform all the relevant error checks.
			{
				EnergyCheck(EnP, EneLog[(PSave-2)], AnsObj, p, MRate, Pshift, EneLog, ParamLog);
				EneLog.erase(EneLog.begin());
				EneLog.push_back(EnP);
				ParamLog.erase(ParamLog.begin());
				ParamLog.push_back(AnsObj->ParamList());				
			}
			if (TermFlag)
			{
				EnIter = EnIter.head(p);
				break;
			}
			if (SampFlag == false)
			{
				pcount++;
				P_fin = P_fin + AnsObj->ParamList();
			}
		}
		P_fin = P_fin / pcount;
		AnsObj->ParamLoad(P_fin);
		return EnIter;
	}
}

#endif