// ------------------------------------------------------------------------------------------------------------
// configuration.hh - contains the declarations of the Config object class, which represents the physical state
//				      of the system at any given point in the Markov chain. The Config class does not contain 
//				      any symmetries of the system - these are governed by the Hilbert object class that is 
//				      used to generate Config objects and supply proposed changes that respect any symmetries.
// ------------------------------------------------------------------------------------------------------------
#ifndef CONFIGURATION_HH
#define CONFIGURATION_HH

#include <iostream>
#include <vector>

using namespace std;

namespace nqsvmc 
{
	// Definition of the Diff struct, included here to be consistent across all configuration types.
	struct Diff
	{
		int num; // - Number of altered sites.
		vector<int> pos; // - Vector of lattice site position indices.
		vector<int> val; // - Vector of local operator value differences.
		int sign; // - Sign change associated with the difference - normally 1, only relevant for fermions.
	};

	class Configuration // - The overarching parent class for all configurations.
	{
	public: // Only virtual public functions included here - definitions are on a subclass basis.
		// Observer functions for extracting parameter values from a generic Config object:
		// - Size will return the size of the lattice this configuration pertains to.
		virtual int Size()const = 0;
		// - SzTot will return the total spin projection - spinless systems return zero.
		virtual int SzTot()const = 0;
		// - Nparticle will return the total number of mobile particles - spin systems return zero.
		virtual int Nparticle()const = 0;
		// - FullCfg will return a vector output representing the entire configuration.
		virtual vector<int> FullCfg()const = 0;
		// Configuration alteration functions:
		// - CfgChange will alter the existing configuration when supplied a Diff struct.
		// - Not to be confused with Diff2Cfg - Diff2Cfg will generate new configurations rather than alter an
		// - existing Config object, and is a method associated with the Hilbert class.
		virtual void CfgChange(Diff& diff) = 0;
	};

	// Forward declaration of implementer subclass:
	class Config : public Configuration
	{
	private:
		Configuration* c_;
	public:
		// Definitions of functions in Config:
		Config(vector<int>& state, char cfgtype);
		// - Constructor for fermionic configurations.
		Config(vector<int>& upstate, vector<int>& dnstate);
		// - Construct a spin or bosonic configuration based on string.
		void SelectType(vector<int>& state, char cfgtype);
		// Observer functions for extracting values from a Config object:
		// - Size will return the number of sites in the configuration.
		int Size() const;
		// - SzTot will return the total spin projection - spinless systems return zero.
		int SzTot()const;
		// - Nparticle will return the total number of mobile particles - spin systems return zero.
		int Nparticle()const;
		// - FullCfg will return a vector output representing the entire configuration.
		vector<int> FullCfg() const;
		// Configuration alteration functions:
		// - CfgChange will alter the existing configuration when supplied a Diff struct.
		void CfgChange(Diff& diff);
	};
}

#include "spinconfig.hh"
#include "boseconfig.hh"
#include "fermconfig.hh"

#endif