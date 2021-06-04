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
}
#endif