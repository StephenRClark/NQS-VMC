// ------------------------------------------------------------------------------------------------------------
// Hilbert.cpp - contains the definitions of the Hilbert class of object. Since the subclasses execute their
//				 methods differently, a general Hilbert object is effectively a wrapper that asks a new 
//				 subclass instance to perform their functions and returns those results.
// ------------------------------------------------------------------------------------------------------------

#ifndef HILBERT_CPP
#define HILBERT_CPP

#include "Hilbert.hh"
#include "Configuration.hh"

namespace nqsvmc
{
	class Spin : public GeneralSpin
	{
	private:
		GeneralSpin* s_;
	public:
		// Generate a populated version of the Bose object when given the necessary arguments.
		Spin(int size, int sztot, double smag, bool fixed)
		{
			if (fixed) // fixed = true will construct the fixed projection Spin subclass.
			{
				s_ = new SpinFZ(size, sztot, smag);
			}
			else
			{
				s_ = new SpinVZ(size, sztot, smag);
			}
		}
		// Observer functions for extracting parameter values will involve calling the new generated subclass instance.
		int SiteDim()const
		{
			return s_->SiteDim();
		}
		int SysSize()const
		{
			return s_->SysSize();
		}
		// Configuration and difference return will depend on configuration restrictions.
		// - RandomCfg will return a configuration of the appropriate subclass.
		Config* RandomCfg() const
		{
			return s_->RandomCfg();
		}
		// - PropMove will return a Diff struct for an appropriate proposed configuration.
		Diff PropMove(Config* Cfg) const
		{
			return s_->PropMove(Cfg);
		}
		// - Diff2Cfg will output new BoseCfg objects related to the input SpinCfg by the supplied Diff structs.
		Config* Diff2Cfg(Config* Cfg, Diff diff)const
		{
			return s_->Diff2Cfg(Cfg, diff);
		}
	};

	class Bose : public GeneralBose
	{
	private:
		GeneralBose* b_;
	public:
		// Generate a populated version of the Bose object when given the necessary arguments.
		Bose(int size, int nbose, int nmax, bool fixed)
		{
			if (fixed) // fixed signifies whether the total number is fixed (true) or not.
			{
				b_ = new BoseFN(size, nbose, nmax);
			}
			else
			{
				b_ = new BoseVN(size, nbose, nmax);
			}
		}
		// Observer functions for extracting parameter values will involve calling the new generated subclass instance.
		int SiteDim()const
		{
			return b_->SiteDim();
		}
		int SysSize()const
		{
			return b_->SysSize();
		}
		// Configuration and difference return will depend on configuration restrictions.
		// - RandomCfg will return a configuration of the appropriate subclass.
		Config* RandomCfg() const
		{
			return b_->RandomCfg();
		}
		// - PropMove will return a Diff struct for an appropriate proposed configuration.
		Diff PropMove(Config* Cfg) const
		{
			return b_->PropMove(Cfg);
		}
		// - Diff2Cfg will output new BoseCfg objects related to the input BoseCfg by the supplied Diff structs.
		Config* Diff2Cfg(Config* Cfg, Diff diff)const
		{
			return b_->Diff2Cfg(Cfg, diff);
		}
	};

	class Ferm : public GeneralFerm
	{
	private:
		GeneralFerm* f_;
	public:
		// Generate a populated version of the Bose object when given the necessary arguments.
		Ferm(int size, int nferm, int sztot, char constraint)
		{
			if (constraint == 'x') // An entry of x will be used to signify no restrictions.
			{
				f_ = new FermVN(size, nferm, sztot);
			}
			else if (constraint == 'n') // An entry of n will signify an overall fixed number.
			{
				f_ = new FermFN(size, nferm, sztot);
			}
			else if (constraint == 's') // An entry of s will signify fixed number and total spin.
			{
				f_ = new FermFZ(size, nferm, sztot);
			}
			else
			{
				cerr << "Invalid constraint identifier. Use 'x' (none), 'n' (fixed number) or 's' (fixed number and spin)." << endl;
				std::abort();
			}
		}
		// Observer functions for extracting parameter values will involve calling the new generated subclass instance.
		int SiteDim()const
		{
			return f_->SiteDim();
		}
		int SysSize()const
		{
			return f_->SysSize();
		}
		// Configuration and difference return will depend on configuration restrictions.
		// - RandomCfg will return a configuration of the appropriate subclass.
		Config* RandomCfg() const
		{
			return f_->RandomCfg();
		}
		// - PropMove will return a Diff struct for an appropriate proposed configuration.
		Diff PropMove(Config* Cfg) const
		{
			return f_->PropMove(Cfg);
		}
		// - Diff2Cfg will output new BoseCfg objects related to the input SpinCfg by the supplied Diff structs.
		Config* Diff2Cfg(Config* Cfg, Diff diff)const
		{
			return f_->Diff2Cfg(Cfg, diff);
		}
	};

	// Definitions of functions in Hilbert:
	// - Constructor for spin Hilbert - third input is a double, for spin magnitude.
	Hilbert::Hilbert(int size, int sztot, double smag, bool constraint)
	{
		h_ = new Spin(size, sztot, smag, constraint);
		typeID = 's';
	}
	// - Constructor for bosonic Hilbert - last input is a boolean flag for fixed / variable number.
	Hilbert::Hilbert(int size, int nbose, int nmax, bool constraint)
	{
		h_ = new Bose(size, nbose, nmax, constraint);
		typeID = 'b';
	}
	// - Constructor for fermionic Hilbert - last input is a character, for different constraint types.
	Hilbert::Hilbert(int size, int nferm, int sztot, char constraint)
	{
		h_ = new Ferm(size, nferm, sztot, constraint);
		typeID = 'f';
	}
	// Observer functions:
	// Type is unique to the implementer subclass, and is used to identify.
	char Hilbert::Type() const
	{
		return typeID;
	}
	int Hilbert::SiteDim() const
	{
		return h_->SiteDim();
	}
	int Hilbert::SysSize() const
	{
		return h_->SysSize();
	}
	// Configuration manipulation functions:
	Config* Hilbert::RandomCfg() const
	{
		return h_->RandomCfg();
	}
	Diff Hilbert::PropMove(Config* Cfg) const
	{
		return h_->PropMove(Cfg);
	}
	Config* Hilbert::Diff2Cfg(Config* Cfg, Diff diff) const
	{
		return h_->Diff2Cfg(Cfg, diff);
	}
}
#endif