// ------------------------------------------------------------------------------------------------------------
// configuration.cpp - contains the definitions of the Config class of object. Since the subclasses execute 
//				       their methods differently, a general Config object is effectively a wrapper that asks a
//					   new subclass instance to perform their functions and returns those results.
// ------------------------------------------------------------------------------------------------------------

#ifndef CONFIGURATION_CPP
#define CONFIGURATION_CPP

#include "Configuration.hh"

namespace nqsvmc
{
	// Definitions for the Config wrapper class.
	Config::Config(vector<int>& state, char cfgtype)
	{
		SelectType(state, cfgtype); // Type selection performed in SelectType.
	}
	// - Constructor for fermionic configurations.
	Config::Config(vector<int>& upstate, vector<int>& dnstate)
	{
		c_ = new FermCfg(upstate, dnstate);
	}
	// - Construct a spin or bosonic configuration based on string.
	void Config::SelectType(vector<int>& state, char cfgtype)
	{
		if (cfgtype == 's')
		{
			c_ = new SpinCfg(state);
		}
		else if (cfgtype == 'b')
		{
			c_ = new BoseCfg(state);
		}
		else
		{
			cerr << "Invalid configuration type identifier - use 's' or 'b'." << endl;
			std::abort();
		}
	}
	// Observer functions for extracting values from a Config object:
	// - Size will return the number of sites in the configuration.
	int Config::Size() const
	{
		return c_->Size();
	}
	// - SzTot will return the total spin projection - spinless systems return zero.
	int Config::SzTot()const
	{
		return c_->SzTot();
	}
	// - Nparticle will return the total number of mobile particles - spin systems return zero.
	int Config::Nparticle()const
	{
		return c_->Nparticle();
	}
	// - FullCfg will return a vector output representing the entire configuration.
	vector<int> Config::FullCfg() const
	{
		return c_->FullCfg();
	}
	// Configuration alteration functions:
	// - CfgChange will alter the existing configuration when supplied a Diff struct.
	void Config::CfgChange(Diff& diff)
	{
		return c_->CfgChange(diff);
	}

	// Definitions of Configuration subclass functions here:

	// <<<<<<----------- BoseCfg functions ------------>>>>>>

	// Constructor functions for SpinCfg subclass:
		// - For no inputs, will assume N = 0, Sz = 0, d = 1;
	SpinCfg::SpinCfg()
	{
		vector<int> v = { 0 };
		N = 0;
		Sz = 0;
		spins = v;
	}
	// - Given an input vector of spins, will create the associated SpinCfg object.
	SpinCfg::SpinCfg(vector<int>& state)
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
	// - SpinCfg::Size will return the number of sites in the configuration.
	int SpinCfg::Size()const
	{
		return N;
	}
	// - SpinCfg::SzTot will return the total Sz projection of the configuration.
	int SpinCfg::SzTot()const
	{
		return Sz;
	}
	// - SpinCfg::Nparticle will return 0 as spins are treated as static objects / features.
	int SpinCfg::Nparticle()const
	{
		return 0;
	}
	// - SpinCfg::FullCfg will return a N x 1 vector of local Sz projections.
	vector<int> SpinCfg::FullCfg()const
	{
		return spins;
	}
	// Configuration alteration functions:
	// - SpinCfg::CfgChange will alter the existing SpinCfg object when supplied a Diff.
	void SpinCfg::CfgChange(Diff& diff)
	{
		// For each altered site, adjust the local spin projection.
		for (int d = 0; d < diff.num; d++)
		{
			spins[diff.pos[d]] += diff.val[d];
			Sz += diff.val[d];
		}
	}

	// <<<<<<----------- BoseCfg functions ------------>>>>>>

	// Constructor functions for BoseCfg subclass:
	// - For no inputs, will assume N = 0, Nb = 0, d = 1.
	BoseCfg::BoseCfg()
	{
		vector<int> v = { 0 };
		N = 0;
		Nb = 0;
		occ = v;
	}
	// - Given an input vector of occupations, will create the associated BoseCfg object.
	BoseCfg::BoseCfg(vector<int>& state)
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
		// - BoseCfg::Size will return the number of sites in the configuration.
	int BoseCfg::Size()const
	{
		return N;
	}
	// - BoseCfg::SzTot will return 0 as boson implementation is currently spinless.
	int BoseCfg::SzTot()const
	{
		return 0;
	}
	// - BoseCfg::Nparticle will return the total number of bosons across the configuration.
	int BoseCfg::Nparticle()const
	{
		return Nb;
	}
	// - BoseCfg::FullCfg will return a N x 1 vector of on-site occupations.
	vector<int> BoseCfg::FullCfg()const
	{
		return occ;
	}
	// Configuration alteration functions:
	// - BoseCfg::CfgChange will alter the existing BoseCfg object when supplied a Diff.
	void BoseCfg::CfgChange(Diff& diff)
	{
		// For each altered site, adjust the on-site occupation.
		for (int d = 0; d < diff.num; d++)
		{
			occ[diff.pos[d]] += diff.val[d];
			Nb += diff.val[d];
		}
		// Going to take it on faith that the Diff values are sensible and don't result in negative values...
	}

	// <<<<<<----------- FermCfg functions ------------>>>>>>

	// Constructor functions for the FermCfg subclass:
		// - For no inputs, will assume N = 0, Nf = 0, d = 1;
	FermCfg::FermCfg()
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
	FermCfg::FermCfg(vector<int>& f_up, vector<int>& f_dn)
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
		ext.resize((size_t)2*N);
		for (int i = 0; i < N; i++)
		{
			N_up += f_up[i];
			N_dn += f_dn[i];
			ext[(size_t)N + i] = f_dn[i]; // Append elements of f_dn to ext to construct full vector.
			// Assuming resize works fine, should not require constant adjustment of storage.
		}
		N_ext = ext;
	}
	// Observer functions for extracting parameter values from a FermCfg object:
	// - FermCfg::Size will return the number of sites in the configuration.
	int FermCfg::Size()const
	{
		return N;
	}
	// - FermCfg::SzTot will return the total Sz projection of the configuration.
	int FermCfg::SzTot()const
	{
		int total = N_up + N_dn;
		return total;
	}
	// - FermCfg::Nparticle will return the total number of fermions in the configuration.
	int FermCfg::Nparticle()const
	{
		int total = N_up - N_dn;
		return total;
	}
	// - FermCfg::FullCfg will return a 2N x 1 vector of occupations, with sites N and greater being down spins
	vector<int> FermCfg::FullCfg()const
	{
		return N_ext;
	}
	// Configuration alteration functions:
	// - FermCfg::CfgChange will alter the existing FermCfg object when supplied a Diff.
	void FermCfg::CfgChange(Diff& diff)
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
}
#endif