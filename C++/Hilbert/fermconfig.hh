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
	public:
		// Constructor functions for the FermCfg subclass:
		// - For no inputs, will assume N = 0, Nf = 0, d = 1;
		FermCfg()
		{
			vector<int> v = { 0 };
			N = 0;
			N_up = 0;
			up = v;
			N_dn = 0;
			dn = v;
			N_ext = v;
		}
		// - Given two vectors of up-spin and down-spin fermions, will create the associated FermCfg object.
		FermCfg(vector<int>& f_up, vector<int>& f_dn)
		{
			N = (int)f_up.size();
			if (N != f_dn.size()) // Check that the input vectors are of the same size.
			{
				cerr << "Up and down fermion occupation vectors are not of equal size." << endl;
				std::abort();
			}
			N_up = 0;
			N_dn = 0;
			up = f_up;
			dn = f_dn;
			vector<int> ext = f_up;
			ext.resize(2 * N);
			for (int i = 0; i < N; i++)
			{
				N_up += f_up[i];
				N_dn += f_dn[i];
				ext[N + i] = f_dn[i]; // Append elements of f_dn to ext to construct full vector.
				// Assuming resize works fine, should not require constant adjustment of storage.
			}
			N_ext = ext;
		}
		// Observer functions for extracting parameter values from a FermCfg object:
		// - Size will return the number of sites in the configuration.
		int Size()const
		{
			return N;
		}
		// - SzTot will return the total Sz projection of the configuration.
		int SzTot()const
		{
			int total = N_up + N_dn;
			return total;
		}
		// - Nparticle will return the total number of fermions in the configuration.
		int Nparticle()const
		{
			int total = N_up - N_dn;
			return total;
		}
		// - FullCfg will return a 2N x 1 vector of occupations, with sites N and greater being down spins
		vector<int> FullCfg()const
		{
			return N_ext;
		}
		// Configuration alteration functions:
		// - CfgChange will alter the existing FermCfg object when supplied a Diff.
		void CfgChange(Diff& diff)
		{
			int site;
			// Diffs will be generated such that the changes to down occupations have 'positions' I = i+N.
			for (int d = 0; d < diff.num; d++)
			{
				site = diff.pos[d] % N; // Down fermion changes listed in positions >= N. This line selects the lattice site.
				N_ext[diff.pos[d]] += diff.val[d]; // Can alter the extended vector straight away.
				if (diff.pos[d] > site) // Down fermion change.
				{
					dn[site] += diff.val[d];
				}
				else // Up fermion change.
				{
					up[site] += diff.val[d];
				}
			}
		}
	private:
		int N; // - Number of sites.
		int N_up; // - Number of up-spin fermions.
		vector<int> up; // - N x 1 ordered vector of sites with up-spin fermions.
		int N_dn; // - Number of down-spin fermions.
		vector<int> dn; // - N x 1 ordered vector of sites with down spin fermions.
		vector<int> N_ext; // - 2N x 1 ordered vector of all sites, with down fermions in sites N+1:2N.
	};	
}

#endif