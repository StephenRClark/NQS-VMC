// ------------------------------------------------------------------------------------------------------------
// Ansatz.cpp - contains the definitions of the Ansatz object class, outlining its interactions with the Config
//			    type of object. The Ansatz class contains all the methods of calculating wavefunction 
//			    properties such as amplitudes and parameter derivatives. Within the Ansatz is a Reference 
//			    object and some number of Modifier objects, and the methods in Ansatz will call the methods of 
//				its constituents.
// ------------------------------------------------------------------------------------------------------------

#ifndef ANSATZ_CPP
#define ANSATZ_CPP

#include <Eigen/Dense>
#include "Ansatz.hh"

using namespace std;
using namespace Eigen;

namespace nqsvmc
{
	// Definitions of functions in Ansatz:
	// Constructor for the Ansatz object type:
	Ansatz::Ansatz(Reference* Ref, vector<Modifier*> Mods, Hilbert* Hlb)
	{ // Inclusion of Hilbert to ensure compatibility of components.
		h_ = Hlb;
		r_ = Ref;
		m_ = Mods;
		Nvar = (int)r_->VarFlag();
		NpTotal = ((int)r_->VarFlag() * r_->Nparam());
		for (int m = 0; m < m_.size(); m++)
		{
			Nvar += (int)m_[m]->VarFlag();
			NpTotal += ((int)m_[m]->VarFlag() * m_[m]->Nparam());
		}
	}

	// Observer functions for extracting information from an Ansatz object:
	// - NumVar returns the number of variational components in the Ansatz.
	int Ansatz::NumVar() const
	{
		return Nvar;
	}
	// - Nparam returns the number of variational parameters in the Ansatz.
	int Ansatz::Nparam() const
	{
		return NpTotal;
	}
	// - Nmod returns the number of Modifiers.
	int Ansatz::Nmod() const
	{
		return (int)m_.size();
	}
	// - ParamList will return a vector of parameters in the ansatz from References and Modifiers.
	VectorXd Ansatz::ParamList() const
	{
		VectorXd P = VectorXd::Zero(NpTotal);
		int p = 0;
		if (r_->VarFlag())
		{
			int npr = r_->Nparam();
			VectorXd dPr = r_->ParamList();
			if (npr == NpTotal)
			{
				P = dPr;
			}
			else
			{
				P.segment(p, npr) = dPr;
			}
			p += npr;
		}
		for (int m = 0; m < m_.size(); m++)
		{
			if (m_[m]->VarFlag())
			{
				int npm = m_[m]->Nparam();
				VectorXd dPm = m_[m]->ParamList();
				if (npm == NpTotal)
				{
					P = dPm;
				}
				else
				{
					P.segment(p, npm) = dPm;
				}
				p += npm;
			}
		}
		return P;
	}

	// Variational management functions:
	// - RndBtcSelect will randomly select parameters within the Ansatz components.
	// -- First version with no specified batch fraction.
	void Ansatz::RndBtcSelect()
	{
		if (r_->VarFlag())
		{
			r_->RndBtcSelect();
		}
		for (int m = 0; m < m_.size(); m++)
		{
			if (m_[m]->VarFlag())
			{
				m_[m]->RndBtcSelect();
			}
		}
		return;
	}
	// -- Second version uses specified batch fraction.
	void Ansatz::RndBtcSelect(double pfrac)
	{
		if (r_->VarFlag())
		{
			r_->RndBtcSelect(pfrac);
		}
		for (int m = 0; m < m_.size(); m++)
		{
			if (m_[m]->VarFlag())
			{
				m_[m]->RndBtcSelect(pfrac);
			}
		}
		return;
	}

	// Wavefunction calculation and update functions:
	// Function here to trim Infs and NaNs before they cause trouble.
	void ParamCheck(VectorXd& P)
	{
		for (int p = 0; p < P.size(); p++)
		{
			if (isnan(P(p)) || isinf(P(p)))
			{
				P(p) = 0;
			}
			if (abs(P(p)) < 1e-30)
			{
				P(p) = 0;
			}
		}
		return;
	}
	// - PsiUpdate will update the variational parameters of the wavefunction when given a vector of changes.
	void Ansatz::PsiUpdate(VectorXd dP)
	{
		ParamCheck(dP);
		int p = 0;
		if (r_->VarFlag())
		{
			VectorXd dPr;
			int npr = r_->Nparam();
			if (npr == NpTotal)
			{
				dPr = dP;
			}
			else
			{
				dPr = dP.segment(0, npr);
			}
			r_->PsiUpdate(dPr);
			p += npr;
		}
		for (int m = 0; m < m_.size(); m++)
		{
			if (m_[m]->VarFlag())
			{
				VectorXd dPm;
				int npm = m_[m]->Nparam();
				if (npm == NpTotal)
				{
					dPm = dP;
				}
				else
				{
					dPm = dP.segment(p, npm);
				}
				m_[m]->PsiUpdate(dPm);
				p += npm;
			}
		}
		return;
	}
	// - ParamLoad will overwrite the variatonal parameters of the wavefunction when given a vector of parameters.
	void Ansatz::ParamLoad(VectorXd P)
	{
		ParamCheck(P);
		if (P.size() != NpTotal)
		{
			cerr << "Parameter vector does not have the correct number of elements." << endl;
			std::abort();
		}
		int p = 0;
		if (r_->VarFlag())
		{
			VectorXd dPr;
			int npr = r_->Nparam();
			if (npr == NpTotal)
			{
				dPr = P;
			}
			else
			{
				dPr = P.segment(0, npr);
			}
			r_->ParamLoad(dPr);
			p += npr;
		}
		for (int m = 0; m < m_.size(); m++)
		{
			if (m_[m]->VarFlag())
			{
				VectorXd dPm;
				int npm = m_[m]->Nparam();
				if (npm == NpTotal)
				{
					dPm = P;
				}
				else
				{
					dPm = P.segment(p, npm);
				}
				m_[m]->ParamLoad(dPm);
				p += npm;
			}
		}
		return;
	}
	// - OptIndLoad will overwrite the variatonal parameters of the wavefunction when given a vector of parameters.
	void Ansatz::OptIndLoad(VectorXd newinds)
	{
		ParamCheck(newinds);
		if (newinds.size() != NpTotal)
		{
			cerr << "Parameter vector does not have the correct number of elements." << endl;
			std::abort();
		}
		int p = 0;
		if (r_->VarFlag())
		{
			VectorXd dPr;
			int npr = r_->Nparam();
			if (npr == NpTotal)
			{
				dPr = newinds;
			}
			else
			{
				dPr = newinds.segment(0, npr);
			}
			r_->OptIndLoad(dPr);
			p += npr;
		}
		for (int m = 0; m < m_.size(); m++)
		{
			if (m_[m]->VarFlag())
			{
				VectorXd dPm;
				int npm = m_[m]->Nparam();
				if (npm == NpTotal)
				{
					dPm = newinds;
				}
				else
				{
					dPm = newinds.segment(p, npm);
				}
				m_[m]->OptIndLoad(dPm);
				p += npm;
			}
		}
		return;
	}
	// - PrepPsi will load any local configuration information used in the wavefunction constituents.
	void Ansatz::PrepPsi(Config* Cfg)
	{
		r_->PrepPsi(Cfg);
		for (int m = 0; m < m_.size(); m++)
		{
			m_[m]->PrepPsi(Cfg);
		}
		return;
	}
	// - PsiCfgUpdate will update any local configuration information after a configuration change.
	void Ansatz::PsiCfgUpdate()
	{
		r_->PsiCfgUpdate();
		for (int m = 0; m < m_.size(); m++)
		{
			m_[m]->PsiCfgUpdate();
		}
		return;
	}
	// Proposed local information changes contained within each Modifier.
	// If required later, will also have PsiCfgUpdateExt which uses an external local information source.
	// - PsiRatio will return the ratio of two amplitudes when supplied the configuration difference and the
	// -- appropriate local information for one configuration loaded into the Ansatz.
	double Ansatz::PsiRatio(Diff diff)
	{
		double Ratio = 1;
		for (int m = 0; m < m_.size(); m++)
		{
			Ratio *= m_[m]->PsiRatio(diff);
		}
		Ratio *= r_->PsiRatio(diff);
		if (isnan(Ratio) || isinf(Ratio))
		{
			Ratio = 0;
		}
		return Ratio;
	}
	// If required later, will also have PsiRatioExt which loads local information to an external object.
	// - LogDeriv will return the logarithmic derivatives of the wavefunction w.r.t. its parameters.
	VectorXd Ansatz::LogDeriv(Config* Cfg)
	{
		VectorXd dP = VectorXd::Zero(NpTotal);
		int p = 0;
		if (r_->VarFlag())
		{
			VectorXd dPr = r_->LogDeriv(Cfg);
			int npr = (int)dPr.size();
			if (npr == NpTotal)
			{
				dP = dPr;
			}
			else
			{
				dP.segment(p, npr) = dPr;
			}
			p += npr;
		}
		for (int m = 0; m < m_.size(); m++)
		{
			if (m_[m]->VarFlag())
			{
				int npm = m_[m]->Nparam();
				VectorXd dPm = m_[m]->LogDeriv(Cfg);
				if (npm == NpTotal)
				{
					dP = dPm;
				}
				else
				{
					dP.segment(p, npm) = dPm;
				}
				p += npm;
			}
		}
		ParamCheck(dP);
		return dP;
	}
}

#endif