// ------------------------------------------------------------------------------------------------------------
// Hilbert.cpp - contains the definition of the Hilbert class of object. Since the subclasses execute their
//				 methods differently, a general Hilbert object is effectively a wrapper that asks a new 
//				 subclass instance to perform their functions and returns those results.
// ------------------------------------------------------------------------------------------------------------

#ifndef HILBERT_CPP
#define HILBERT_CPP

#include "Hilbert.hh"
#include "Configuration.hh"

namespace nqsvmc
{
	// Definition of generic Spin wrapper class.
	class Spin : public GeneralSpin
	{
	private:
		GeneralSpin* s_;
	public:
		// Generate a populated version of the Bose object when given the necessary arguments.
		Spin(int size, int sztot, double smag, bool fixed)
		{
			if (fixed) // fixed = true will construct the fixed projection Spin subclass.
			{
				s_ = new SpinFZ(size, sztot, smag);
			}
			else
			{
				s_ = new SpinVZ(size, sztot, smag);
			}
		}
		// Observer functions for extracting parameter values will involve calling the new generated subclass instance.
		int SiteDim()const
		{
			return s_->SiteDim();
		}
		int SysSize()const
		{
			return s_->SysSize();
		}
		// Configuration and difference return will depend on configuration restrictions.
		// - RandomCfg will return a configuration of the appropriate subclass.
		Config* RandomCfg() const
		{
			return s_->RandomCfg();
		}
		// - PropMove will return a Diff struct for an appropriate proposed configuration.
		Diff PropMove(Config* Cfg) const
		{
			return s_->PropMove(Cfg);
		}
		// - Diff2Cfg will output new BoseCfg objects related to the input SpinCfg by the supplied Diff structs.
		Config* Diff2Cfg(Config* Cfg, Diff diff)const
		{
			return s_->Diff2Cfg(Cfg, diff);
		}
	};

	// Definition of generic Bose wrapper class.
	class Bose : public GeneralBose
	{
	private:
		GeneralBose* b_;
	public:
		// Generate a populated version of the Bose object when given the necessary arguments.
		Bose(int size, int nbose, int nmax, bool fixed)
		{
			if (fixed) // fixed signifies whether the total number is fixed (true) or not.
			{
				b_ = new BoseFN(size, nbose, nmax);
			}
			else
			{
				b_ = new BoseVN(size, nbose, nmax);
			}
		}
		// Observer functions for extracting parameter values will involve calling the new generated subclass instance.
		int SiteDim()const
		{
			return b_->SiteDim();
		}
		int SysSize()const
		{
			return b_->SysSize();
		}
		// Configuration and difference return will depend on configuration restrictions.
		// - RandomCfg will return a configuration of the appropriate subclass.
		Config* RandomCfg() const
		{
			return b_->RandomCfg();
		}
		// - PropMove will return a Diff struct for an appropriate proposed configuration.
		Diff PropMove(Config* Cfg) const
		{
			return b_->PropMove(Cfg);
		}
		// - Diff2Cfg will output new BoseCfg objects related to the input BoseCfg by the supplied Diff structs.
		Config* Diff2Cfg(Config* Cfg, Diff diff)const
		{
			return b_->Diff2Cfg(Cfg, diff);
		}
	};
	// Definition of generic Ferm wrapper class.
	class Ferm : public GeneralFerm
	{
	private:
		GeneralFerm* f_;
	public:
		// Generate a populated version of the Bose object when given the necessary arguments.
		Ferm(int size, int nferm, int sztot, char constraint)
		{
			if (constraint == 'x') // An entry of x will be used to signify no restrictions.
			{
				f_ = new FermVN(size, nferm, sztot);
			}
			else if (constraint == 'n') // An entry of n will signify an overall fixed number.
			{
				f_ = new FermFN(size, nferm, sztot);
			}
			else if (constraint == 's') // An entry of s will signify fixed number and total spin.
			{
				f_ = new FermFZ(size, nferm, sztot);
			}
			else
			{
				cerr << "Invalid constraint identifier. Use 'x' (none), 'n' (fixed number) or 's' (fixed number and spin)." << endl;
				std::abort();
			}
		}
		// Observer functions for extracting parameter values will involve calling the new generated subclass instance.
		int SiteDim()const
		{
			return f_->SiteDim();
		}
		int SysSize()const
		{
			return f_->SysSize();
		}
		// Configuration and difference return will depend on configuration restrictions.
		// - RandomCfg will return a configuration of the appropriate subclass.
		Config* RandomCfg() const
		{
			return f_->RandomCfg();
		}
		// - PropMove will return a Diff struct for an appropriate proposed configuration.
		Diff PropMove(Config* Cfg) const
		{
			return f_->PropMove(Cfg);
		}
		// - Diff2Cfg will output new Cfg objects related to the input Cfg by the supplied Diff structs.
		Config* Diff2Cfg(Config* Cfg, Diff diff)const
		{
			return f_->Diff2Cfg(Cfg, diff);
		}
	};

	// Definitions of functions in Hilbert wrapper class:
	// - Constructor for spin Hilbert - third input is a double, for spin magnitude.
	Hilbert::Hilbert(int size, int sztot, double smag, bool constraint)
	{
		h_ = new Spin(size, sztot, smag, constraint);
		typeID = 's';
	}
	// - Constructor for bosonic Hilbert - last input is a boolean flag for fixed / variable number.
	Hilbert::Hilbert(int size, int nbose, int nmax, bool constraint)
	{
		h_ = new Bose(size, nbose, nmax, constraint);
		typeID = 'b';
	}
	// - Constructor for fermionic Hilbert - last input is a character, for different constraint types.
	Hilbert::Hilbert(int size, int nferm, int sztot, char constraint)
	{
		h_ = new Ferm(size, nferm, sztot, constraint);
		typeID = 'f';
	}
	// Observer functions:
	// Type is unique to the implementer subclass, and is used to identify.
	char Hilbert::Type() const
	{
		return typeID;
	}
	int Hilbert::SiteDim() const
	{
		return h_->SiteDim();
	}
	int Hilbert::SysSize() const
	{
		return h_->SysSize();
	}
	// Configuration manipulation functions:
	Config* Hilbert::RandomCfg() const
	{
		return h_->RandomCfg();
	}
	Diff Hilbert::PropMove(Config* Cfg) const
	{
		return h_->PropMove(Cfg);
	}
	Config* Hilbert::Diff2Cfg(Config* Cfg, Diff diff) const
	{
		return h_->Diff2Cfg(Cfg, diff);
	}

	// Definitions of Hilbert subclass functions here:

	// <<<<<<------------ GeneralSpin functions ------------>>>>>>

	// - Diff2Cfg will output new SpinCfg objects related to the input SpinCfg by the supplied Diff structs.
	Config* GeneralSpin::Diff2Cfg(Config* Cfg, Diff diff) const
	{
		vector<int> state(Cfg->FullCfg());
		for (int d = 0; d < diff.num; d++)
		{
			state[diff.pos[d]] += diff.val[d];
		}
		char type = 's';
		Config* NewConf = new Config(state, type);
		return NewConf;
	}

	// <<<<<<------------ Spin subclass functions ------------>>>>>>

	// --> SpinFZ
	// Constructor for the fixed total Spin Hilbert subclass.
	SpinFZ::SpinFZ(int size, int sztot, double smag)
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
	// - SpinFZ::SiteDim returns the single site local dimension of the Hilbert space.
	int SpinFZ::SiteDim()const
	{
		return int(2 * S + 1);
	}
	// - SpinFZ::SysSize returns the number of sites in the system.
	int SpinFZ::SysSize()const
	{
		return N;
	}
	// SpinFZ::RandomCfg will generate a configuration with the total spin specified in Hilbert.
	Config* SpinFZ::RandomCfg() const
	{
		vector<int> state{ 0 };
		state.resize(N);
		random_device rd;
		std::default_random_engine rnd(rd());
		uniform_int_distribution<int> sitegen(0, N - 1);
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
		Config* Conf = new Config(state, type);
		return Conf;
	}
	// - SpinFZ::PropMove will lower one spin and raise another.
	Diff SpinFZ::PropMove(Config* Cfg) const
	{
		Diff diff;
		diff.num = 2;
		diff.sign = 1;
		diff.pos.resize(2); diff.val.resize(2);
		diff.val[0] = -2; diff.val[1] = 2;
		vector<int> state(Cfg->FullCfg());
		int N_bot = 0;
		int dN_bot = 0;
		for (int i = 0; i < N; i++)
		{
			N_bot += (int)(state[N] == (int)(-2 * S));
		}
		random_device rd;
		std::default_random_engine rnd(rd());
		uniform_int_distribution<int> sitegen(0, N - 1);
		diff.pos[0] = sitegen(rnd);
		diff.pos[1] = sitegen(rnd);
		while (state[diff.pos[0]] > (-2 * S)) // Ensure a valid site is selected first.
		{
			diff.pos[0] = sitegen(rnd);
		}
		int S_0 = state[diff.pos[0]];
		dN_bot += (int)(S_0 == (int)(-2 * S + 1));
		while (state[diff.pos[1]] >= (2 * S)) // Ensure destination is not a site already at max projection.
		{
			diff.pos[1] = sitegen(rnd);
		}
		int S_P = state[diff.pos[1]];
		dN_bot -= (int)(S_P == (int)(-2 * S));
		double N_botP = (double)N_bot + dN_bot;
		diff.tfac = ((double)N - N_botP) / ((double)N - (double)N_bot);
		return diff;
	}

	// --> SpinVZ
	// Constructor for the variable total Spin Hilbert subclass.
	SpinVZ::SpinVZ(int size, int sztot, double smag)
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
	// - SpinVZ::SiteDim returns the single site local dimension of the Hilbert space.
	int SpinVZ::SiteDim()const
	{
		return int(2 * S + 1);
	}
	// - SpinVZ::SysSize returns the number of sites in the system.
	int SpinVZ::SysSize()const
	{
		return N;
	}
	// - SpinVZ::RandomCfg will generate a configuration with random spin projections.
	Config* SpinVZ::RandomCfg() const
	{
		vector<int> state{ 0 };
		uniform_int_distribution<int> spingen(0, int(2 * S));
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
		Config* Conf = new Config(state, type);
		return Conf;
	}
	// - SpinVZ::PropMove will alter the spin of one site randomly.
	Diff SpinVZ::PropMove(Config* Cfg) const
	{
		Diff diff;
		diff.num = 1;
		diff.sign = 1;
		diff.pos.resize(1);
		diff.val.resize(1);
		diff.tfac = 1.0;
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
		int value = int(newval - oldval);
		diff.pos[0] = site;
		diff.val[0] = value;
		return diff;
	}

	// <<<<<<------------ GeneralBose functions ------------>>>>>>
	
	// - Diff2Cfg will output new BoseCfg objects related to the input BoseCfg by the supplied Diff structs.
	Config* GeneralBose::Diff2Cfg(Config* Cfg, Diff diff)const
	{
		vector<int> state(Cfg->FullCfg());
		for (int d = 0; d < diff.num; d++)
		{
			state[diff.pos[d]] += diff.val[d];
		}
		char type = 'b';
		Config* NewConf = new Config(state, type);
		return NewConf;
	};

	// <<<<<<------------ Bose subclass functions ------------>>>>>>

	// --> BoseFN
	// Constructor for the fixed number Bose Hilbert subclass.
	BoseFN::BoseFN(int size, int nbose, int nmax)
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
	// - BoseFN::SiteDim returns the single site local dimension of the Hilbert space.
	int BoseFN::SiteDim() const
	{
		return Nmax + 1;
	}
	// - BoseFN::SysSize returns the number of sites in the system.
	int BoseFN::SysSize() const
	{
		return N;
	}
	// - BoseFN::RandomCfg will propose a configuration according to the number of bosons specified.
	Config* BoseFN::RandomCfg() const
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
		Config* Conf = new Config(state, type);
		return Conf;
	}
	// - BoseFN::PropMove will move a boson from one site to another eligible one.
	Diff BoseFN::PropMove(Config* Cfg) const
	{
		Diff diff;
		diff.num = 2;
		diff.sign = 1;
		diff.pos.resize(2); diff.val.resize(2);
		diff.val[0] = -1; diff.val[1] = 1;
		vector<int> state(Cfg->FullCfg());
		int N_occ = 0;
		int dN_occ = 0;
		for (int i = 0; i < N; i++)
		{
			N_occ += (int)(state[i] > 0);
		}
		uniform_int_distribution<int> sitegen(0, N - 1);
		random_device rd;
		std::default_random_engine rnd(rd());
		diff.pos[0] = sitegen(rnd);
		diff.pos[1] = sitegen(rnd);
		while (state[diff.pos[0]] == 0) // Ensure a populated site is selected first.
		{
			diff.pos[0] = sitegen(rnd);
		}
		int N_0 = state[diff.pos[0]];
		dN_occ -= (int)(N_0 == 1);
		// Ensure destination is not a site already at max capacity or starting site.
		while ((state[diff.pos[1]] >= Nmax) || (diff.pos[0] == diff.pos[1]))
		{
			diff.pos[1] = sitegen(rnd);
		}
		int N_P = state[diff.pos[1]];
		dN_occ += (int)(N_P == 0);
		double N_occP = (double)N_occ + dN_occ;
		diff.tfac = N_occP / ((double)N_occ);
		return diff;
	}

	// --> BoseVN
	// Constructor for the variable number Bose Hilbert subclass.
	BoseVN::BoseVN(int size, int nbose, int nmax)
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
	// - BoseVN::SiteDim returns the single site local dimension of the Hilbert space.
	int BoseVN::SiteDim() const
	{
		return Nmax + 1;
	}
	// - BoseVN::SysSize returns the number of sites in the system.
	int BoseVN::SysSize() const
	{
		return N;
	}
	// BoseVN::RandomCfg will generate an entirely random starting configuration.
	Config* BoseVN::RandomCfg() const
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
		Config* Conf = new Config(state, type);
		return Conf;
	};
	// BoseVN::PropMove will change the occupation of one of the sites.
	Diff BoseVN::PropMove(Config* Cfg) const
	{
		Diff diff;
		diff.num = 1;
		diff.sign = 1;
		diff.pos.resize(1);
		diff.val.resize(1);
		diff.tfac = 1.0;
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

	// <<<<<<------------ GeneralFerm functions ------------>>>>>>

	// - Diff2Cfg will output new BoseCfg objects related to the input BoseCfg by the supplied Diff structs.
	Config* GeneralFerm::Diff2Cfg(Config* Cfg, Diff diff) const
	{
		vector<int> state(Cfg->FullCfg());
		size_t N = state.size();
		for (int d = 0; d < diff.num; d++)
		{
			state[diff.pos[d]] += diff.val[d];
		}
		vector<int> upstate;
		vector<int> dnstate;
		upstate.resize(N);
		dnstate.resize(N);
		for (size_t i = 0; i < N; i++)
		{
			upstate[i] = state[i];
			dnstate[i] = state[i + N];
		}
		Config* NewConf = new Config(upstate, dnstate);
		return NewConf;
	}

	// <<<<<<------------ Ferm subclass functions ------------>>>>>>

	// --> FermFZ
	// Constructor for the variable number Ferm Hilbert subclass.
	FermFZ::FermFZ(int size, int nferm, int sztot)
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
	// - FermFZ::SiteDim returns the single site local dimension of the Hilbert space.
	int FermFZ::SiteDim() const
	{
		return 4;
	}
	// - FermFZ::SysSize returns the number of sites in the system.
	int FermFZ::SysSize() const
	{
		return N;
	}
	// - FermFZ::RandomCfg will generate a configuration with the specified start population.
	Config* FermFZ::RandomCfg() const
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
	// - FermFZ::PropMove will move a fermion from one site to another, spin specific.
	Diff FermFZ::PropMove(Config* Cfg) const
	{
		vector<int> state(Cfg->FullCfg());
		uniform_int_distribution<int> longdist(0, 2 * N - 1);
		uniform_int_distribution<int> shrtdist(0, N - 1);
		random_device rd;
		std::default_random_engine rnd(rd());
		Diff diff;
		diff.num = 2; // Two site occupations are altered.
		diff.val.resize(2);
		diff.val[0] = -1; diff.val[1] = 1;
		diff.pos[0] = longdist(rnd);
		diff.tfac = 1.0;
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

	// --> FermFN
	// Constructor for the variable number Ferm Hilbert subclass.
	FermFN::FermFN(int size, int nferm, int sztot)
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
	// - FermFN::SiteDim returns the single site local dimension of the Hilbert space.
	int FermFN::SiteDim() const
	{
		return 4;
	}
	// - FermFN::SysSize returns the number of sites in the system.
	int FermFN::SysSize() const
	{
		return N;
	}
	// - FermFN::RandomCfg will generate a configuration with the specified start population.
	Config* FermFN::RandomCfg() const
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
	// - FermFN::PropMove will move a fermion from one site to another, spin agnostic.
	Diff FermFN::PropMove(Config* Cfg) const
	{
		vector<int> state(Cfg->FullCfg());
		uniform_int_distribution<int> longdist(0, 2 * N - 1);
		random_device rd;
		std::default_random_engine rnd(rd());
		Diff diff;
		diff.num = 2; // Two site occupations are altered.
		diff.val.resize(2);
		diff.val[0] = -1; diff.val[1] = 1;
		diff.pos.resize(2);
		diff.pos[0] = longdist(rnd);
		diff.pos[1] = longdist(rnd);
		diff.tfac = 1.0;
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

	// --> FermVN
	// Constructor for the variable number Ferm Hilbert subclass.
	FermVN::FermVN(int size, int nferm, int sztot)
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
	// - FermVN::SiteDim returns the single site local dimension of the Hilbert space.
	int FermVN::SiteDim() const
	{
		return 4;
	}
	// - FermFN::SysSize returns the number of sites in the system.
	int FermVN::SysSize() const
	{
		return N;
	}
	// - FermVN::RandomCfg will generate a configuration with the specified start population.
	Config* FermVN::RandomCfg() const
	{
		vector<int> upstate{ 0 }, dnstate{ 0 };
		uniform_int_distribution<int> shrtdist(0, N - 1);
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
	// - FermVN::PropMove will alter occupation of one random site.
	Diff FermVN::PropMove(Config* Cfg) const
	{
		vector<int> state(Cfg->FullCfg());
		uniform_int_distribution<int> longdist(0, 2 * N - 1);
		random_device rd;
		std::default_random_engine rnd(rd());
		int site = longdist(rnd);
		Diff diff;
		diff.num = 1;
		diff.tfac = 1.0;
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
}
#endif