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
	private:
		int N; // - Number of sites.
		int Nb; // - Total number of bosons.
		vector<int> occ; // - Vector of on-site occupation numbers.
	public:
		// Constructor functions for BoseCfg subclass:
		// - For no inputs, will assume N = 0, Nb = 0, d = 1.
		BoseCfg();
		// - Given an input vector of occupations, will create the associated BoseCfg object.
		BoseCfg(vector<int>& state);
		// Observer functions for extracting parameter values from a BoseCfg object:
		// - Size will return the number of sites in the configuration.
		int Size()const;
		// - SzTot will return the total spin projection. Boson implementation is currently spinless.
		int SzTot()const;
		// - Nparticle will return the total number of bosons across the configuration.
		int Nparticle()const;
		// - FullCfg will return a N x 1 vector of on-site occupations.
		vector<int> FullCfg()const;
		// Configuration alteration functions:
		// - CfgChange will alter the existing BoseCfg object when supplied a Diff.
		void CfgChange(Diff& diff);	
	};
}

#endif