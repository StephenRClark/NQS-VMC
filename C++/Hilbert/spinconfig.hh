// ------------------------------------------------------------------------------------------------------------
// spinconfig.hh - contains the definitions of the SpinCfg subclass of Config objects.
// ------------------------------------------------------------------------------------------------------------

#ifndef SPINCONFIG_HH
#define SPINCONFIG_HH

#include <vector>
#include "Configuration.hh"

using namespace std;

namespace nqsvmc
{
	class SpinCfg : public Configuration
	{
	private:
		int N; // - Number of sites.
		int Sz; // - Total spin z-projection.
		vector<int> spins; // - Vector of local spin z-projections.
	public:
		// Constructor functions for SpinCfg subclass:
		// - For no inputs, will assume N = 0, Sz = 0, d = 1;
		SpinCfg();
		// - Given an input vector of spins, will create the associated SpinCfg object.
		SpinCfg(vector<int>& state);
		// Observer functions for extracting parameter values from a SpinCfg object:
		// - Size will return the number of sites in the configuration.
		int Size()const;
		// - SzTot will return the total Sz projection of the configuration.
		int SzTot()const;
		// - Nparticle will return the number of particles in the configuration. Spins are treated as static objects / features.
		int Nparticle()const;
		// - FullCfg will return a N x 1 vector of local Sz projections.
		vector<int> FullCfg()const;
		// Configuration alteration functions:
		// - CfgChange will alter the existing SpinCfg object when supplied a Diff.
		void CfgChange(Diff& diff);
	};
}

#endif