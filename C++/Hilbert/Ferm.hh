// ------------------------------------------------------------------------------------------------------------
// Ferm.hh - contains the definitions of the Ferm class of Hilbert object, detailing its interactions with the 
//			 FermCfg type of Config object. The Hilbert class of object is intended for setting rules on what
//			 configurations are permitted, and generating Config objects subject to those rules.
// ------------------------------------------------------------------------------------------------------------

#ifndef FERM_HH
#define FERM_HH

#include <vector>
#include "Hilbert.hh"
#include "Configuration.hh"

using namespace std;

namespace nqsvmc 
{
	class GeneralFerm : public GeneralHilbert
	{
	public:
		virtual int SiteDim() const = 0;
		virtual int SysSize() const = 0;
		virtual Config* RandomCfg() const = 0;
		virtual Diff PropMove(Config* Cfg) const = 0;
		virtual Config* Diff2Cfg(Config* Cfg, Diff diff) const = 0;
	};
	

	class FermVN : public GeneralFerm // Subclass for variable total number and spin.
	{
	private:
		int N; // Number of sites.
		int N_up; // Number of up fermions.
		int N_dn; // Number of down fermions.
	public:
		// Constructor for the variable number Ferm Hilbert subclass.
		FermVN(int size, int nferm, int sztot)
		{
			if (abs(sztot) > nferm)
			{
				cerr << "Desired spin projection is larger than permitted by fermion number." << endl;
				std::abort();
			}
			if (nferm > 2 * size)
			{
				cerr << "Number of fermions is larger than permitted by system size." << endl;
				std::abort();
			}
			if ((nferm % 2) != (sztot % 2))
			{
				cerr << "Desired fermion count and spin projection are incompatible." << endl;
				std::abort();
			}
			if (size <= 0)
			{
				cerr << "Invalid system size." << endl;
				std::abort();
			}
			N = nferm;
			N_up = (nferm + sztot) / 2;
			N_dn = (nferm - sztot) / 2;
		}
		// Observer functions:
		int SiteDim() const
		{
			return 4;
		}
		int SysSize() const
		{
			return N;
		}
		// Configuration manipulation functions:
		// - RandomCfg will generate a configuration with the populations specified in Hilbert.
		Config* RandomCfg() const
		{
			vector<int> upstate{ 0 }, dnstate{ 0 };
			uniform_int_distribution<int> shrtdist(0, N-1);
			upstate.resize(N);
			dnstate.resize(N);
			random_device rd;
			std::default_random_engine rnd(rd());
			int site;
			for (int i = 0; i < N_up; i++)
			{
				site = shrtdist(rnd);
				while (upstate[site] == 1) // Ensure target site is not already populated.
				{
					site = shrtdist(rnd);
				}
				upstate[site] = 1;
			}
			for (int j = 0; j < N_dn; j++)
			{
				site = shrtdist(rnd);
				while (dnstate[site] == 1) // Ensure target site is not already populated.
				{
					site = shrtdist(rnd);
				}
				dnstate[site] = 1;
			}
			Config* Conf = new Config(upstate, dnstate);
			return Conf;
		}
		// PropMove will alter occupation of one random site.
		Diff PropMove(Config* Cfg) const
		{
			vector<int> state(Cfg->FullCfg());
			uniform_int_distribution<int> longdist(0, 2*N - 1);
			random_device rd;
			std::default_random_engine rnd(rd());
			int site = longdist(rnd);
			Diff diff;
			diff.num = 1;
			diff.val.resize(1);
			diff.pos.resize(1);
			diff.val[0] = 1 - 2 * state[site];
			diff.pos[0] = site;
			state[site] += 1;
			state[site] = state[site] % 2;
			int nswap = 0; // Counter for operator swaps invoked with this move.
			for (int i = 0; i < site; i++)
			{
				nswap += state[i];
			}
			diff.sign = 1 - 2 * (nswap % 2);
			return diff;
		}
		Config* Diff2Cfg(Config* Cfg, Diff diff) const
		{
			vector<int> state(Cfg->FullCfg());
			for (int d = 0; d < diff.num; d++)
			{
				state[diff.pos[d]] += diff.val[d];
			}
			vector<int> upstate;
			vector<int> dnstate;
			upstate.resize(N);
			dnstate.resize(N);
			for (int i = 0; i < N; i++)
			{
				upstate[i] = state[i];
				dnstate[i] = state[i + N];
			}
			Config* NewConf = new Config(upstate, dnstate);
			return NewConf;
		}
	};

	class FermFN : public GeneralFerm // Subclass for fixed total number and variable spin.
	{
	private:
		int N; // Number of sites.
		int N_up; // Number of up fermions.
		int N_dn; // Number of down fermions.
	public:
		// Constructor for the variable number Ferm Hilbert subclass.
		FermFN(int size, int nferm, int sztot)
		{
			if (abs(sztot) > nferm)
			{
				cerr << "Desired spin projection is larger than permitted by fermion number." << endl;
				std::abort();
			}
			if (nferm > 2 * size)
			{
				cerr << "Number of fermions is larger than permitted by system size." << endl;
				std::abort();
			}
			if ((nferm % 2) != (sztot % 2))
			{
				cerr << "Desired fermion count and spin projection are incompatible." << endl;
				std::abort();
			}
			if (size <= 0)
			{
				cerr << "Invalid system size." << endl;
				std::abort();
			}
			N = nferm;
			N_up = (nferm + sztot) / 2;
			N_dn = (nferm - sztot) / 2;
		}
		// Observer functions:
		int SiteDim() const
		{
			return 4;
		}
		int SysSize() const
		{
			return N;
		}
		// Configuration manipulation functions:
		// - RandomCfg will generate a configuration with the populations specified in Hilbert.
		Config* RandomCfg() const
		{
			vector<int> upstate{ 0 }, dnstate{ 0 };
			uniform_int_distribution<int> shrtdist(0, N - 1);
			random_device rd;
			std::default_random_engine rnd(rd());
			upstate.resize(N);
			dnstate.resize(N);
			int site;
			for (int i = 0; i < N_up; i++)
			{
				site = shrtdist(rnd);
				while (upstate[site] == 1) // Ensure target site is not already populated.
				{
					site = shrtdist(rnd);
				}
				upstate[site] = 1;
			}
			for (int j = 0; j < N_dn; j++)
			{
				site = shrtdist(rnd);
				while (dnstate[site] == 1) // Ensure target site is not already populated.
				{
					site = shrtdist(rnd);
				}
				dnstate[site] = 1;
			}
			Config* Conf = new Config(upstate, dnstate);
			return Conf;
		}
		// PropMove will move a fermion from one site to another, spin agnostic.
		Diff PropMove(Config* Cfg) const
		{
			vector<int> state(Cfg->FullCfg());
			uniform_int_distribution<int> longdist(0, 2*N - 1);
			random_device rd;
			std::default_random_engine rnd(rd());
			Diff diff;
			diff.num = 2; // Two site occupations are altered.
			diff.val.resize(2);
			diff.val[0] = -1; diff.val[1] = 1;
			diff.pos.resize(2);
			diff.pos[0] = longdist(rnd);
			diff.pos[1] = longdist(rnd);
			while (state[diff.pos[0]] == 0) // Ensure starting location is occupied.
			{
				diff.pos[0] = longdist(rnd);
			}
			while (state[diff.pos[1]] == 1) // Ensure destination is not occupied.
			{
				diff.pos[1] = longdist(rnd);
			}
			int nswap = 0; // Counter for operator swaps invoked with this move.
			int site1 = min(diff.pos[0], diff.pos[1]);
			int site2 = max(diff.pos[0], diff.pos[1]);
			for (int i = site1; i < site2; i++)
			{
				nswap += state[i];
			}
			diff.sign = 1 - 2 * (nswap % 2);
			return diff;
		}
		Config* Diff2Cfg(Config* Cfg, Diff diff) const
		{
			vector<int> state(Cfg->FullCfg());
			for (int d = 0; d < diff.num; d++)
			{
				state[diff.pos[d]] += diff.val[d];
			}
			vector<int> upstate;
			vector<int> dnstate;
			upstate.resize(N);
			dnstate.resize(N);
			for (int i = 0; i < N; i++)
			{
				upstate[i] = state[i];
				dnstate[i] = state[i + N];
			}
			Config* NewConf = new Config(upstate, dnstate);
			return NewConf;
		}
	};

	class FermFZ : public GeneralFerm // Subclass for fixed total number and spin.
	{
	private:
		int N; // Number of sites.
		int N_up; // Number of up fermions.
		int N_dn; // Number of down fermions.
	public:
		// Constructor for the variable number Ferm Hilbert subclass.
		FermFZ(int size, int nferm, int sztot)
		{
			if (abs(sztot) > nferm)
			{
				cerr << "Desired spin projection is larger than permitted by fermion number." << endl;
				std::abort();
			}
			if (nferm > 2 * size)
			{
				cerr << "Number of fermions is larger than permitted by system size." << endl;
				std::abort();
			}
			if ((nferm % 2) != (sztot % 2))
			{
				cerr << "Desired fermion count and spin projection are incompatible." << endl;
				std::abort();
			}
			if (size <= 0)
			{
				cerr << "Invalid system size." << endl;
				std::abort();
			}
			N = nferm;
			N_up = (nferm + sztot) / 2;
			N_dn = (nferm - sztot) / 2;
		}
		// Observer functions:
		int SiteDim() const
		{
			return 4;
		}
		int SysSize() const
		{
			return N;
		}
		// Configuration manipulation functions:
		// - RandomCfg will generate a configuration with the populations specified in Hilbert.
		Config* RandomCfg() const
		{
			vector<int> upstate{ 0 }, dnstate{ 0 };
			uniform_int_distribution<int> shrtdist(0, N - 1);
			random_device rd;
			std::default_random_engine rnd(rd());
			upstate.resize(N);
			dnstate.resize(N);
			int site;
			for (int i = 0; i < N_up; i++)
			{
				site = shrtdist(rnd);
				while (upstate[site] == 1) // Ensure target site is not already populated.
				{
					site = shrtdist(rnd);
				}
				upstate[site] = 1;
			}
			for (int j = 0; j < N_dn; j++)
			{
				site = shrtdist(rnd);
				while (dnstate[site] == 1) // Ensure target site is not already populated.
				{
					site = shrtdist(rnd);
				}
				dnstate[site] = 1;
			}
			Config* Conf = new Config(upstate, dnstate);
			return Conf;
		}
		// - PropMove will move a fermion from one site to another, spin specific.
		Diff PropMove(Config* Cfg) const
		{
			vector<int> state(Cfg->FullCfg());
			uniform_int_distribution<int> longdist(0, 2*N - 1);
			uniform_int_distribution<int> shrtdist(0, N - 1);
			random_device rd;
			std::default_random_engine rnd(rd());
			Diff diff;
			diff.num = 2; // Two site occupations are altered.
			diff.val.resize(2);
			diff.val[0] = -1; diff.val[1] = 1;
			diff.pos[0] = longdist(rnd);
			while (state[diff.pos[0]] == 0) // Ensure starting location is occupied.
			{
				diff.pos[0] = longdist(rnd);
			}
			int smod = 0;
			if (diff.pos[0] > N)
			{
				smod = N;
			}
			diff.pos[1] = shrtdist(rnd) + smod;
			while (state[diff.pos[1]] == 1)
			{
				diff.pos[1] = shrtdist(rnd) + smod;
			}
			int nswap = 0; // Counter for operator swaps invoked with this move.
			int site1 = min(diff.pos[0], diff.pos[1]);
			int site2 = max(diff.pos[0], diff.pos[1]);
			for (int i = site1; i < site2; i++) // Sum over occupied sites between the two positions.
			{
				nswap += state[i];
			}
			diff.sign = 1 - 2 * (nswap % 2);
			return diff;
		}
		Config* Diff2Cfg(Config* Cfg, Diff diff) const
		{
			vector<int> state(Cfg->FullCfg());
			for (int d = 0; d < diff.num; d++)
			{
				state[diff.pos[d]] += diff.val[d];
			}
			vector<int> upstate;
			vector<int> dnstate;
			upstate.resize(N);
			dnstate.resize(N);
			for (int i = 0; i < N; i++)
			{
				upstate[i] = state[i];
				dnstate[i] = state[i + N];
			}
			Config* NewConf = new Config(upstate, dnstate);
			return NewConf;
		}
	};
}
#endif