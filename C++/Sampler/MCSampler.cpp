// ------------------------------------------------------------------------------------------------------------
// MCSampler.cpp - contains the implementation of the MCSampler object class, which executes Markov Chain Monte
//				   Carlo sampling with a given Ansatz and Hamiltonian.
// ------------------------------------------------------------------------------------------------------------

#ifndef MCSAMPLER_CPP
#define MCSAMPLER_CPP

#include <Eigen/Dense>
#include "MCSampler.hh"

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	// Constructor functions:
		// - Basic usage, only assign necessary fields and use defaults for the rest.
	MCSampler::MCSampler(Hilbert* hlb, Hamiltonian* hmt)
	{
		hb_ = hlb;
		hm_ = hmt;
	}
	// - Advanced usage, assign all parameters at initialisation.
	MCSampler::MCSampler(Hilbert* hlb, Hamiltonian* hmt, int Nsample, int Nequilibrate, int Nbl)
	{
		hb_ = hlb;
		hm_ = hmt;
		Nsamp = Nsample;
		Nequil = Nequilibrate;
		Nblock = Nbl;
	}
	// Observer functions:
	// - NumSamp returns the number of sample steps.
	int MCSampler::NumSamp() const
	{
		return Nsamp;
	}
	// - NumBlock returns the number of steps between samples.
	int MCSampler::NumBlock() const
	{
		return Nblock;
	}

	// Parameter assignment functions:
	// - SetNsamp is used to modify the number of samples.
	void MCSampler::SetNsamp(int Nnew)
	{
		if (Nnew <= 0)
		{
			cerr << "Number of samples cannot be zero or negative.";
			std::abort();
		}
		Nsamp = Nnew;
	}
	// - SetNblock is used to modify the number of steps per sample.
	void MCSampler::SetNblock(int Nnew)
	{
		if (Nnew <= 0)
		{
			cerr << "Number of sample steps cannot be zero or negative.";
			std::abort();
		}
		Nblock = Nnew;
	}
	// - SetNequil is used to modify the number of equilibration steps before sampling.
	void MCSampler::SetNequil(int Nnew)
	{
		if (Nnew < 0)
		{
			cerr << "Number of equilibration steps cannot be negative.";
			std::abort();
		}
		Nequil = Nnew;
	}
	// - SetHilbert will replace the Hilbert pointer with a new one.
	void MCSampler::SetHilbert(Hilbert* hnew)
	{
		hb_ = hnew;
		return;
	}
	// - SetHamiltonian will replace the Hamiltonian pointer with a new one.
	void MCSampler::SetHamiltonian(Hamiltonian* hnew)
	{
		hm_ = hnew;
		return;
	}
	// Sampling functions:
	// - MCMCSample - sample the Ansatz object and output local energy, derivatives and move rate.
	void MCSampler::MCMCSample(Ansatz* AnsObj, double& EnAvg, VectorXd& dP, vector<Operator*> opptrs_,
		vector<MatrixXd>& vals, double& MRate) const
	{
		// Initialise parameters and storage.
		if (opptrs_.size() != vals.size())
		{
			cerr << "Size of value storage (" << vals.size() << ") does not match number of operators ("
				<< opptrs_.size() << ")." << endl;
			std::abort();
		}
		int NP = AnsObj->Nparam();
		int N_op = (int)opptrs_.size();
		EnAvg = 0;
		MRate = 0;
		dP.setZero();
		Diff diff;
		double Ratio = 0;
		double TProb = 0;
		double EnLoc = 0;
		double R = 0;
		VectorXd dLogp = VectorXd::Zero(NP);
		random_device rd;
		default_random_engine rndgen(rd());
		uniform_real_distribution<double> randdist(0, 1);
		Config* cfg = hb_->RandomCfg();
		AnsObj->PrepPsi(cfg);
		for (int o = 0; o < N_op; o++)
		{
			MatrixXd locval = opptrs_[o]->GraphSample(cfg, EnLoc, dLogp, AnsObj);
			vals[o] = locval.setZero();
		}
		for (int q = 0; q < Nequil; q++)
		{
			diff = hb_->PropMove(cfg);
			Ratio = AnsObj->PsiRatio(diff);
			TProb = pow(Ratio, 2)*diff.tfac;
			R = randdist(rndgen);
			if (R < TProb)
			{
				cfg->CfgChange(diff);
				AnsObj->PsiCfgUpdate();
			}
		}
		for (int ns = 0; ns < Nsamp; ns++)
		{
			for (int nb = 0; nb < Nblock; nb++)
			{
				diff = hb_->PropMove(cfg);
				Ratio = AnsObj->PsiRatio(diff);
				TProb = pow(Ratio, 2)*diff.tfac;
				R = randdist(rndgen);
				if (R < TProb)
				{
					cfg->CfgChange(diff);
					AnsObj->PsiCfgUpdate();
					MRate += (double)(1 / ((double)Nsamp * (double)Nblock));
				}
			}
			EnLoc = hm_->EnergySample(cfg, AnsObj);
			if (isnan(EnLoc) || isinf(EnLoc))
			{
				EnLoc = 0; // Trying some forward error correction to avoid NaN's.
			}
			EnAvg += (EnLoc / (double)Nsamp);
			dLogp = AnsObj->LogDeriv(cfg);
			dP += (dLogp / (double)Nsamp);
			for (int o = 0; o < N_op; o++)
			{
				vals[o] += (opptrs_[o]->GraphSample(cfg, EnLoc, dLogp, AnsObj) / (double)Nsamp);
			}
		}
		return;
	}
	// - EvalSample - sample the Ansatz object and output local energy.
	void MCSampler::EvalSample(Ansatz* AnsObj, double& EnAvg, vector<Operator*> opptrs_, vector<MatrixXd>& vals) const
	{
		// Initialise parameters and storage.
		if (opptrs_.size() != vals.size())
		{
			cerr << "Size of value storage (" << vals.size() << ") does not match number of operators ("
				<< opptrs_.size() << ")." << endl;
			std::abort();
		}
		int NP = AnsObj->Nparam();
		int N_op = (int)opptrs_.size();
		EnAvg = 0;
		Diff diff;
		double Ratio = 0;
		double TProb = 0;
		double EnLoc = 0;
		double R = 0;
		VectorXd dLogp = VectorXd::Zero(NP);
		random_device rd;
		default_random_engine rndgen(rd());
		uniform_real_distribution<double> randdist(0, 1);
		Config* cfg = hb_->RandomCfg();
		AnsObj->PrepPsi(cfg);
		for (int o = 0; o < N_op; o++)
		{
			MatrixXd locval = opptrs_[o]->GraphSample(cfg, EnLoc, dLogp, AnsObj);
			vals[o] = locval.setZero();
		}
		for (int q = 0; q < Nequil; q++)
		{
			diff = hb_->PropMove(cfg);
			Ratio = AnsObj->PsiRatio(diff);
			TProb = pow(Ratio, 2)*diff.tfac;
			R = randdist(rndgen);
			if (R < TProb)
			{
				cfg->CfgChange(diff);
				AnsObj->PsiCfgUpdate();
			}
		}
		for (int ns = 0; ns < Nsamp; ns++)
		{
			for (int nb = 0; nb < Nblock; nb++)
			{
				diff = hb_->PropMove(cfg);
				Ratio = AnsObj->PsiRatio(diff);
				TProb = pow(Ratio, 2)*diff.tfac;
				R = randdist(rndgen);
				if (R < TProb)
				{
					cfg->CfgChange(diff);
					AnsObj->PsiCfgUpdate();
				}
			}
			EnLoc = hm_->EnergySample(cfg, AnsObj);
			EnAvg += (EnLoc / (double)Nsamp);
			dLogp = AnsObj->LogDeriv(cfg);
			for (int o = 0; o < N_op; o++)
			{
				vals[o] += (opptrs_[o]->GraphSample(cfg, EnLoc, dLogp, AnsObj) / (double)Nsamp);
			}
		}
		return;
	}
}

#endif