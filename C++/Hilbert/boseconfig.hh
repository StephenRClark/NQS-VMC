// ------------------------------------------------------------------------------------------------------------
// boseconfig.hh - contains the definitions of the BoseCfg subclass of the Config object type.
// ------------------------------------------------------------------------------------------------------------

#ifndef BOSECONFIG_HH
#define BOSECONFIG_HH

#include <vector>
#include "Configuration.hh"

using namespace std;

namespace nqsvmc 
{
	class BoseCfg : public Configuration
	{
	public:
		// Constructor functions for BoseCfg subclass:
		// - For no inputs, will assume N = 0, Nb = 0, d = 1.
		BoseCfg()
		{
			vector<int> v = { 0 };
			N = 0;
			Nb = 0;
			occ = v;
		}
		// - Given an input vector of occupations, will create the associated BoseCfg object.
		BoseCfg(vector<int>& state)
		{
			N = (int)state.size();
			// occ.resize(N);
			Nb = 0;
			for (int i = 0; i < N; i++)
			{
				if (state[i] < 0) // Check there are no negative numbers in input state.
				{
					cerr << "Input state has negative occupation numbers." << endl;
					std::abort();
				}
				Nb += state[i];
				// occ[i] = state[i];
			}
			occ = state;
		}
		// Observer functions for extracting parameter values from a BoseCfg object:
		// - Size will return the number of sites in the configuration.
		int Size()const 
		{
			return N;
		}
		// - SzTot will return the total spin projection. Boson implementation is currently spinless.
		int SzTot()const
		{
			return 0;
		}
		// - Nparticle will return the total number of bosons across the configuration.
		int Nparticle()const
		{
			return Nb;
		}
		// - FullCfg will return a N x 1 vector of on-site occupations.
		vector<int> FullCfg()const
		{
			return occ;
		}
		// Configuration alteration functions:
		// - CfgChange will alter the existing BoseCfg object when supplied a Diff.
		void CfgChange(Diff& diff) 
		{
			// For each altered site, adjust the on-site occupation.
			for (int d = 0; d < diff.num; d++)
			{
				occ[diff.pos[d]] += diff.val[d];
				Nb += diff.val[d];
			}
			// Going to take it on faith that the Diff values are sensible and don't result in negative values...
		}
	private:
		int N; // - Number of sites.
		int Nb; // - Total number of bosons.
		vector<int> occ; // - Vector of on-site occupation numbers.
	};
}

#endif