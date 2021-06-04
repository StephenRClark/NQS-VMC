// ------------------------------------------------------------------------------------------------------------
// Hilbert.hh - contains the declarations of the Hilbert object class, outlining its interactions with the 
//				Config type of object. The Hilbert class is intended for setting rules on what configurations
//				are permitted, and generating Config objects subject to those rules.
// ------------------------------------------------------------------------------------------------------------
#ifndef HILBERT_HH
#define HILBERT_HH

#include <vector>
#include "Configuration.hh"

using namespace std;

namespace nqsvmc 
{
	class GeneralHilbert // The abstract class that prototypes the general capabilities of a Hilbert object.
	{
	public:
		// Observer functions for extracting parameter values from a generic Hilbert object:
		// - SiteDim will return the on-site dimension of the Hilbert space.
		virtual int SiteDim()const = 0;
		// - SysSize will return the number of sites of the Hilbert space.
		virtual int SysSize()const = 0;
		// Configuration manipulation functions:
		// - RandomCfg will return a configuration of the appropriate subclass.
		virtual Config* RandomCfg() const = 0;
		// - PropMove will return a change in the configuration as a Diff struct.
		// -- Supply the Diff to a CfgChange function to alter the configuration.
		virtual Diff PropMove(Config* Cfg) const = 0;
		// - Diff2Cfg will return configurations when supplied a starting configuration and Diff structs.
		// -- Note that this function is intended to generate new configurations, not alter the existing one.
		// -- Each configuration subclass will be given a separate function that changes the existing Cfg when
		// -- given a singular Diff struct.
		virtual Config* Diff2Cfg(Config* Cfg, Diff diff)const = 0;
	};

	// Forward declaration of implementer subclass:
	class Hilbert : public GeneralHilbert
	{
	private:
		GeneralHilbert* h_;
		char typeID; // Character used as an identifier to check compatibility.
	public:
		// Constructors for the Hilbert class of object.
		// - Constructor for spin Hilbert - third input is a double, for spin magnitude.
		Hilbert(int size, int sztot, double smag, bool constraint);
		// - Constructor for bosonic Hilbert - last input is a boolean flag for fixed / variable number.
		Hilbert(int size, int nbose, int nmax, bool constraint);
		// - Constructor for fermionic Hilbert - last input is a character, for different constraint types.
		Hilbert(int size, int nferm, int sztot, char constraint);
		// Observer functions:
		// Type is unique to the implementer subclass, and is used to identify.
		char Type() const;
		int SiteDim() const;
		int SysSize() const;
		// Configuration manipulation functions:
		Config* RandomCfg() const;
		Diff PropMove(Config* Cfg) const;
		Config* Diff2Cfg(Config* Cfg, Diff diff) const;
	};
}

#include "Spin.hh"
#include "Bose.hh"
#include "Ferm.hh"

#endif