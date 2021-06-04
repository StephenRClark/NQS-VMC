// ------------------------------------------------------------------------------------------------------------
// Spin.hh - contains the definitions of the Spin class of Hilbert object, detailing its interactions with the 
//			 SpinCfg type of Config object. The Hilbert class of object is intended for setting rules on what
//			 configurations are permitted, and generating Config objects subject to those rules.
// ------------------------------------------------------------------------------------------------------------
#ifndef SPIN_HH
#define SPIN_HH

#include <vector>
#include <random>
#include "Hilbert.hh"
#include "Configuration.hh"

using namespace std;

namespace nqsvmc
{
	class GeneralSpin : public GeneralHilbert
	{
	public:
		virtual int SiteDim() const = 0;
		virtual int SysSize() const = 0;
		virtual Config* RandomCfg() const = 0;
		virtual Diff PropMove(Config* Cfg) const = 0;
		virtual Config* Diff2Cfg(Config* Cfg, Diff diff) const = 0;
	};
	class SpinFZ : public GeneralSpin // Subclass for fixed total spin proposed moves.
	{
	private:
		int N; // Number of sites.
		int Sz; // Total z-projection of spin.
		double S; // Magnitude of a single spin.
	public:
		// Constructor for the fixed total Spin Hilbert subclass.
		SpinFZ(int size, int sztot, double smag)
		{
			N = size;
			if (abs(sztot) > (N * smag))
			{
				cerr << "Proposed total spin is larger than permitted." << endl;
				std::abort();
			}
			if (floor(2 * S) != 2 * S)
			{
				cerr << "Require integer or half-integer spin magnitude." << endl;
				std::abort();
			}
			if (smag <= 0)
			{
				cerr << "Invalid spin magnitude." << endl;
				std::abort();
			}
			if (size <= 0)
			{
				cerr << "Invalid system size." << endl;
				std::abort();
			}
			Sz = sztot;
			S = smag;
		}
		// Observer functions:
		int SiteDim()const
		{
			return int(2 * S + 1);
		}
		int SysSize()const
		{
			return N;
		}
		// Configuration manipulation functions:
		// RandomCfg will generate a configuration with the total spin specified in Hilbert.	
		Config* RandomCfg() const
		{
			vector<int> state{ 0 };
			state.resize(N);
			random_device rd;
			std::default_random_engine rnd(rd());
			uniform_int_distribution<int> sitegen(0, N-1);
			for (int i = 0; i < N; i++) // Start with all down state.
			{
				state[i] = int(-2 * S);
			}
			for (int i = 0; i < S * N + Sz; i++)
			{
				int site = sitegen(rnd);
				while (state[site] >= 2 * S) // Ensure target site is not at max Sz value.
				{
					int site = sitegen(rnd);
				}
				state[site] += 2;
			}
			char type = 's';
			Config* Conf = new Config(state,type);
			return Conf;
		}
		// - PropMove will lower one spin and raise another.
		Diff PropMove(Config* Cfg) const
		{
			Diff diff;
			diff.num = 2;
			diff.sign = 1;
			diff.pos.resize(2); diff.val.resize(2);
			diff.val[0] = -2; diff.val[1] = 2;
			vector<int> state(Cfg->FullCfg());
			random_device rd;
			std::default_random_engine rnd(rd());
			uniform_int_distribution<int> sitegen(0, N - 1);
			diff.pos[0] = sitegen(rnd);
			diff.pos[1] = sitegen(rnd);
			while (state[diff.pos[0]] > -2 * S) // Ensure a valid site is selected first.
			{
				diff.pos[0] = sitegen(rnd);
			}
			while (state[diff.pos[1]] >= 2 * S) // Ensure destination is not a site already at max projection.
			{
				diff.pos[1] = sitegen(rnd);
			}			
			return diff;
		}
		Config* Diff2Cfg(Config* Cfg, Diff diff)const
		{
			vector<int> state(Cfg->FullCfg());
			for (int d = 0; d < diff.num; d++)
			{
				state[diff.pos[d]] += diff.val[d];
			}
			char type = 's';
			Config* NewConf = new Config(state,type);
			return NewConf;
		}
	};

	class SpinVZ : public GeneralSpin // Subclass for variable total number proposed moves.
	{
	private:
		int N; // Number of sites.
		int Sz; // Total z-projection of spin.
		double S; // Magnitude of a single spin.
	public:
		// Constructor for the variable total Spin Hilbert subclass.
		SpinVZ(int size, int sztot, double smag)
		{
			N = size;
			if (abs(sztot) > (N * smag))
			{
				cerr << "Proposed total spin is larger than permitted." << endl;
				std::abort();
			}
			if (floor(2 * smag) != 2 * smag)
			{
				cerr << "Require integer or half-integer spin magnitude." << endl;
				std::abort();
			}
			if (smag <= 0)
			{
				cerr << "Invalid spin magnitude." << endl;
				std::abort();
			}
			if (size <= 0)
			{
				cerr << "Invalid system size." << endl;
				std::abort();
			}
			Sz = sztot;
			S = smag;
		}
		// Observer functions:
		int SiteDim()const
		{
			return int(2 * S + 1);
		}
		int SysSize()const
		{
			return N;
		}
		// Configuration manipulation functions:	
		// - RandomCfg will generate a configuration with the total spin specified in Hilbert.
		Config* RandomCfg() const
		{
			vector<int> state{ 0 };
			uniform_int_distribution<int> spingen(0,int(2*S));
			int sztotal = 0;
			state.resize(N);
			random_device rd;
			std::default_random_engine rnd(rd());
			for (int i = 0; i < N; i++) // Randomly assign values to each site.
			{
				state[i] = spingen(rnd) - int(2 * S);
				sztotal += state[i];
			}
			char type = 's';
			Config* Conf = new Config(state,type);
			return Conf;
		}
		// - PropMove will alter one spin randomly.
		Diff PropMove(Config* Cfg) const
		{
			Diff diff;
			diff.num = 1;
			diff.sign = 1;
			diff.pos.resize(1);
			diff.val.resize(1);			
			vector<int> state(Cfg->FullCfg());
			random_device rd;
			std::default_random_engine rnd(rd());
			uniform_int_distribution<int> sitegen(0, N - 1);
			int site = sitegen(rnd);
			int oldval = state[site];
			double newval;
			if (S == (1 / 2))
			{
				newval = -oldval;
			}
			else
			{
				uniform_int_distribution<int> spingen(0, int(2 * S));
				newval = spingen(rnd) - 2 * S;
				while (newval == oldval)
				{
					newval = spingen(rnd) - 2 * S; // Ensure the new value is different from the old.
				}
			}			
			int value = int (newval - oldval);
			diff.pos[0] = site;
			diff.val[0] = value;
			return diff;
		}
		Config* Diff2Cfg(Config* Cfg, Diff diff)const
		{
			vector<int> state(Cfg->FullCfg());
			for (int d = 0; d < diff.num; d++)
			{
				state[diff.pos[d]] += diff.val[d];
			}
			char type = 's';
			Config* NewConf = new Config(state,type);
			return NewConf;
		}
	};
}

#endif