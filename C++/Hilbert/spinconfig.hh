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
	public:
		// Constructor functions for SpinCfg subclass:
		// - For no inputs, will assume N = 0, Sz = 0, d = 1;
		SpinCfg()
		{
			vector<int> v = { 0 };
			N = 0;
			Sz = 0;
			spins = v;
		}
		// - Given an input vector of spins, will create the associated SpinCfg object.
		SpinCfg(vector<int>& state)
		{
			N = (int)state.size();
			Sz = 0;
			for (int i = 0; i < N; i++)
			{
				Sz += state[i];
			}
			spins = state;
		}
		// Observer functions for extracting parameter values from a SpinCfg object:
		// - Size will return the number of sites in the configuration.
		int Size()const
		{
			return N;
		}
		// - SzTot will return the total Sz projection of the configuration.
		int SzTot()const
		{
			return Sz;
		}
		// - Nparticle will return the number of particles in the configuration. Spins are treated as static objects / features.
		int Nparticle()const
		{
			return 0;
		}
		// - FullCfg will return a N x 1 vector of local Sz projections.
		vector<int> FullCfg()const
		{
			return spins;
		}
		// Configuration alteration functions:
		// - CfgChange will alter the existing SpinCfg object when supplied a Diff.
		void CfgChange(Diff& diff)
		{
			// For each altered site, adjust the local spin projection.
			for (int d = 0; d < diff.num; d++)
			{
				spins[diff.pos[d]] += diff.val[d];
				Sz += diff.val[d];
			}
		}
	private:
		int N; // - Number of sites.
		int Sz; // - Total spin z-projection.
		vector<int> spins; // - Vector of local spin z-projections.
	};
}

#endif