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
		virtual int SiteDim() const = 0;
		virtual int SysSize() const = 0;
		virtual Config* RandomCfg() const = 0;
		virtual Diff PropMove(Config* Cfg) const = 0;
		virtual Config* Diff2Cfg(Config* Cfg, Diff diff) const = 0;
	};
	
	class BoseFN : public GeneralBose // Subclass for fixed total number proposed moves.
	{
	private:
		int N; // Number of sites.
		int Nb; // Number of bosons.
		int Nmax; // Maximum occupation of a single site.
	public:
		// Constructor for the fixed number Bose Hilbert subclass.
		BoseFN(int size, int nbose, int nmax)
		{
			N = size;
			if (nbose > (N * nmax))
			{
				cerr << "Proposed number of bosons is larger than permitted." << endl;
				std::abort();
			}
			if (nmax <= 0)
			{
				cerr << "Invalid maximum occupation." << endl;
				std::abort();
			}
			if (size <= 0)
			{
				cerr << "Invalid system size." << endl;
				std::abort();
			}
			Nb = nbose;
			Nmax = nmax;
		}
		// Observer functions:
		int SiteDim()const
		{
			return Nmax + 1;
		}
		int SysSize()const
		{
			return N;
		}
		// Configuration manipulation functions:
		// -RandomCfg will generate a configuration with the number of bosons specified in Hilbert.	
		Config* RandomCfg() const
		{
			vector<int> state{ 0 };
			uniform_int_distribution<int> sitegen(0, N - 1);
			int site;
			state.resize(N);
			random_device rd;
			std::default_random_engine rnd(rd());
			for (int i = 0; i < Nb; i++)
			{
				site = sitegen(rnd);
				while (state[site] >= Nmax)
				{
					site = sitegen(rnd);
				}
				state[site] += 1;
			}
			char type = 'b';
			Config* Conf = new Config(state,type);
			return Conf;
		}
		// PropMove will move a boson from one site to another valid site.
		Diff PropMove(Config* Cfg) const
		{
			vector<int> positions{ 0, 0 };
			vector<int> values{ -1, 1 };
			vector<int> state(Cfg->FullCfg());
			uniform_int_distribution<int> sitegen(0, N - 1);
			random_device rd;
			std::default_random_engine rnd(rd());
			positions[0] = sitegen(rnd);
			positions[1] = sitegen(rnd);
			while (state[positions[0]] == 0) // Ensure a populated site is selected first.
			{
				positions[0] = sitegen(rnd);
			}
			// Ensure destination is not a site already at max capacity or starting site.
			while ((state[positions[1]] >= Nmax) || (positions[0] == positions[1])) 
			{
				positions[1] = sitegen(rnd);
			}
			Diff diff{2,positions,values,1};
			return diff;
		}
		Config* Diff2Cfg(Config* Cfg, Diff diff)const
		{
			vector<int> state(Cfg->FullCfg());
			for (int d = 0; d < diff.num; d++)
			{
				state[diff.pos[d]] += diff.val[d];
			}
			char type = 'b';
			Config* NewConf = new Config(state,type);
			return NewConf;
		}
	};

	class BoseVN : public GeneralBose // Subclass for variable total number proposed moves.
	{
	private:
		int N; // Number of sites.
		int Nb; // Number of bosons.
		int Nmax; // Maximum occupation of a single site.
	public:
		// Constructor for the variable number Bose Hilbert subclass.
		BoseVN(int size, int nbose, int nmax)
		{
			N = size;
			if (nbose > (N * nmax))
			{
				cerr << "Proposed number of bosons is larger than permitted." << endl;
				std::abort();
			}

			if (size <= 0)
			{
				cerr << "Invalid system size." << endl;
				std::abort();
			}
			Nb = nbose;
			Nmax = nmax;
		}
		// Observer functions:
		int SiteDim()const
		{
			return Nmax + 1;
		}
		int SysSize()const
		{
			return N;
		}
		// Configuration manipulation functions:
		// RandomCfg will generate an entirely random starting configuration.
		Config* RandomCfg() const
		{
			vector<int> state{ 0 };
			int nbtot = 0;
			uniform_int_distribution<int> numgen(0, Nmax);
			random_device rd;
			std::default_random_engine rnd(rd());
			state.resize(N);
			for (int i = 0; i < N; i++)
			{
				state[i] = numgen(rnd);
				nbtot += state[i];
			}
			char type = 'b';
			Config* Conf = new Config(state,type);
			return Conf;
		};
		// PropMove will change the occupation of one of the sites.
		Diff PropMove(Config* Cfg) const
		{
			Diff diff;
			diff.num = 1;
			diff.sign = 1;
			diff.pos.resize(1);
			diff.val.resize(1);
			uniform_int_distribution<int> sitegen(0, N - 1);
			uniform_int_distribution<int> numgen(0, Nmax);
			random_device rd;
			std::default_random_engine rnd(rd());
			int site = sitegen(rnd);
			vector<int> state(Cfg->FullCfg());
			int num = numgen(rnd);
			while (state[site] == num) // Ensure proposed occupation is different from existing.
			{
				num = numgen(rnd);
			}			
			diff.pos[0] = site;
			diff.val[0] = num - state[site];			
			return diff;
		};
		Config* Diff2Cfg(Config* Cfg, Diff diff)const
		{
			vector<int> state(Cfg->FullCfg());
			for (int d = 0; d < diff.num; d++)
			{
				state[diff.pos[d]] += diff.val[d];
			}
			char type = 'b';
			Config* NewConf = new Config(state,type);
			return NewConf;
		}
	};
}

#endif