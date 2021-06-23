// ------------------------------------------------------------------------------------------------------------
// Spin.hh - contains the definitions of the Spin class of Hilbert object, detailing its interactions with the 
//			 SpinCfg type of Config object. The Hilbert class of object is intended for setting rules on what
//			 configurations are permitted, and generating Config objects subject to those rules.
// ------------------------------------------------------------------------------------------------------------
#ifndef SPIN_HH
#define SPIN_HH

#include <vector>
#include <random>
#include "Hilbert.hh"
#include "Configuration.hh"

using namespace std;

namespace nqsvmc
{
	class GeneralSpin : public GeneralHilbert
	{
	public:
		// - SiteDim returns the single site local dimension of the Hilbert space.
		virtual int SiteDim() const = 0;
		// - SysSize returns the number of sites in the system.
		virtual int SysSize() const = 0;
		// - RandomCfg will generate a starting random configuration. Subclass specific.
		virtual Config* RandomCfg() const = 0;
		// - PropMove will generate a new configuration from an existing one. Subclass specific.
		virtual Diff PropMove(Config* Cfg) const = 0;
		// - Diff2Cfg will output new SpinCfg objects related to the input SpinCfg by the supplied Diff structs.
		Config* Diff2Cfg(Config* Cfg, Diff diff) const;
	};

	class SpinFZ : public GeneralSpin // Subclass for fixed total spin proposed moves.
	{
	private:
		int N; // Number of sites.
		int Sz; // Total z-projection of spin.
		double S; // Magnitude of a single spin.
	public:
		// Constructor for the fixed total Spin Hilbert subclass.
		SpinFZ(int size, int sztot, double smag);
		// Observer functions need to be overwritten here to access subclass properties.
		// - SiteDim returns the single site local dimension of the Hilbert space.
		int SiteDim() const;
		// - SysSize returns the number of sites in the system.
		int SysSize() const;
		// Configuration manipulation functions:
		// SpinFZ::RandomCfg will generate a configuration with the total spin specified in Hilbert.	
		Config* RandomCfg() const;
		// - PropMove will lower one spin and raise another.
		Diff PropMove(Config* Cfg) const;
		// - Diff2Cfg inherited from GeneralSpin.
	};

	class SpinVZ : public GeneralSpin // Subclass for variable total number proposed moves.
	{
	private:
		int N; // Number of sites.
		int Sz; // Total z-projection of spin.
		double S; // Magnitude of a single spin.
	public:
		// Constructor for the variable total Spin Hilbert subclass.
		SpinVZ(int size, int sztot, double smag);
		// Observer functions need to be overwritten here to access subclass properties.
		// - SiteDim returns the single site local dimension of the Hilbert space.
		int SiteDim() const;
		// - SysSize returns the number of sites in the system.
		int SysSize() const;
		// Configuration manipulation functions:	
		// - SpinVZ::RandomCfg will generate a configuration with random spin projections.
		Config* RandomCfg() const;
		// - SpinVZ::PropMove will alter the spin of one site randomly.
		Diff PropMove(Config* Cfg) const;
		// - Diff2Cfg inherited from GeneralSpin.
	};
}

#endif