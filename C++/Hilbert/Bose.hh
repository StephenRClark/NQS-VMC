// ------------------------------------------------------------------------------------------------------------
// Bose.hh - contains the definitions of the Bose class of Hilbert object, detailing its interactions with the 
//			 BoseCfg type of Config object. The Hilbert class of object is intended for setting rules on what
//			 configurations are permitted, and generating Config objects subject to those rules.
// ------------------------------------------------------------------------------------------------------------
#ifndef BOSE_HH
#define BOSE_HH

#include <vector>
#include <random>
#include "Hilbert.hh"
#include "Configuration.hh"

using namespace std;

namespace nqsvmc
{
	class GeneralBose : public GeneralHilbert
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
		// - Diff2Cfg will output new BoseCfg objects related to the input BoseCfg by the supplied Diff structs.
		Config* Diff2Cfg(Config* Cfg, Diff diff) const;
	};
	
	class BoseFN : public GeneralBose // Subclass for fixed total number proposed moves.
	{
	private:
		int N; // Number of sites.
		int Nb; // Number of bosons.
		int Nmax; // Maximum occupation of a single site.
	public:
		// Constructor for the fixed number Bose Hilbert subclass.
		BoseFN(int size, int nbose, int nmax);
		// Observer functions need to be overwritten here to access subclass properties.
		// - SiteDim returns the single site local dimension of the Hilbert space.
		int SiteDim() const;
		// - SysSize returns the number of sites in the system.
		int SysSize() const;
		// Configuration manipulation functions:
		// - BoseFN::RandomCfg will propose a configuration according to the number of bosons specified.
		Config* RandomCfg() const; 
		// - BoseFN::PropMove will move a boson from one site to another eligible one.
		Diff PropMove(Config* Cfg) const;
		// - Diff2Cfg inherited from GeneralBose.
	};

	class BoseVN : public GeneralBose // Subclass for variable total number proposed moves.
	{
	private:
		int N; // Number of sites.
		int Nb; // Number of bosons.
		int Nmax; // Maximum occupation of a single site.
	public:
		// Constructor for the variable number Bose Hilbert subclass.
		BoseVN(int size, int nbose, int nmax);
		// Observer functions need to be overwritten here to access subclass properties.
		// - SiteDim returns the single site local dimension of the Hilbert space.
		int SiteDim() const;
		// - SysSize returns the number of sites in the system.
		int SysSize() const;
		// Configuration manipulation functions:
		// BoseVN::RandomCfg will generate an entirely random starting configuration.
		Config* RandomCfg() const;
		// BoseVN::PropMove will change the occupation of one of the sites.
		Diff PropMove(Config* Cfg) const;
		// - Diff2Cfg inherited from GeneralBose.
	};
}

#endif