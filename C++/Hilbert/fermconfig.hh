// ------------------------------------------------------------------------------------------------------------
// fermconfig.hh - contains the definitions of the FermCfg subclass of Config objects.
// ------------------------------------------------------------------------------------------------------------

#ifndef FERMCONFIG_HH
#define FERMCONFIG_HH

#include <vector>
#include "Configuration.hh"

using namespace std;

namespace nqsvmc
{
	class FermCfg : public Configuration
	{
	private:
		int N; // - Number of sites.
		int N_up; // - Number of up-spin fermions.
		vector<int> up; // - N x 1 ordered vector of sites with up-spin fermions.
		int N_dn; // - Number of down-spin fermions.
		vector<int> dn; // - N x 1 ordered vector of sites with down spin fermions.
		vector<int> N_ext; // - 2N x 1 ordered vector of all sites, with down fermions in sites N+1:2N.
	public:
		// Constructor functions for the FermCfg subclass:
		// - For no inputs, will assume N = 0, Nf = 0, d = 1;
		FermCfg();
		// - Given two vectors of up-spin and down-spin fermions, will create the associated FermCfg object.
		FermCfg(vector<int>& f_up, vector<int>& f_dn);
		// Observer functions for extracting parameter values from a FermCfg object:
		// - Size will return the number of sites in the configuration.
		int Size()const;
		// - SzTot will return the total Sz projection of the configuration.
		int SzTot()const;
		// - Nparticle will return the total number of fermions in the configuration.
		int Nparticle()const;
		// - FullCfg will return a 2N x 1 vector of occupations, with sites N and greater being down spins
		vector<int> FullCfg()const;
		// Configuration alteration functions:
		// - CfgChange will alter the existing FermCfg object when supplied a Diff.
		void CfgChange(Diff& diff);	
	};	
}

#endif