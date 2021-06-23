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
		// - SiteDim returns the single site local dimension of the Hilbert space.
		virtual int SiteDim() const = 0;
		// - SysSize returns the number of sites in the system.
		virtual int SysSize() const = 0;
		// - RandomCfg will generate a starting random configuration. Subclass specific.
		virtual Config* RandomCfg() const = 0;
		// - PropMove will generate a new configuration from an existing one. Subclass specific.
		virtual Diff PropMove(Config* Cfg) const = 0;
		// - Diff2Cfg will output new SpinCfg objects related to the input SpinCfg by the supplied Diff structs.
		Config* Diff2Cfg(Config* Cfg, Diff diff) const;
	};
	

	class FermVN : public GeneralFerm // Subclass for variable total number and spin.
	{
	private:
		int N; // Number of sites.
		int N_up; // Number of up fermions.
		int N_dn; // Number of down fermions.
	public:
		// Constructor for the variable number Ferm Hilbert subclass.
		FermVN(int size, int nferm, int sztot);
		// Observer functions need to be overwritten here to access subclass properties.
		// - SiteDim returns the single site local dimension of the Hilbert space.
		int SiteDim() const;
		// - SysSize returns the number of sites in the system.
		int SysSize() const;
		// Configuration manipulation functions:
		// - FermVN::RandomCfg will generate a configuration with the specified start population.
		Config* RandomCfg() const;
		// - FermVN::PropMove will alter occupation of one random site.
		Diff PropMove(Config* Cfg) const;
		// - Diff2Cfg inherited from GeneralFerm.		
	};

	class FermFN : public GeneralFerm // Subclass for fixed total number and variable spin.
	{
	private:
		int N; // Number of sites.
		int N_up; // Number of up fermions.
		int N_dn; // Number of down fermions.
	public:
		// Constructor for the variable number Ferm Hilbert subclass.
		FermFN(int size, int nferm, int sztot);
		// Observer functions need to be overwritten here to access subclass properties.
		// - SiteDim returns the single site local dimension of the Hilbert space.
		int SiteDim() const;
		// - SysSize returns the number of sites in the system.
		int SysSize() const;
		// Configuration manipulation functions:
		// - FermFN::RandomCfg will generate a configuration with the specified start population.
		Config* RandomCfg() const;
		// - FermFN::PropMove will move a fermion from one site to another, spin agnostic.
		Diff PropMove(Config* Cfg) const;
		// - Diff2Cfg inherited from GeneralFerm.
	};

	class FermFZ : public GeneralFerm // Subclass for fixed total number and spin.
	{
	private:
		int N; // Number of sites.
		int N_up; // Number of up fermions.
		int N_dn; // Number of down fermions.
	public:
		// Constructor for the variable number Ferm Hilbert subclass.
		FermFZ(int size, int nferm, int sztot);
		// Observer functions need to be overwritten here to access subclass properties.
		// - SiteDim returns the single site local dimension of the Hilbert space.
		int SiteDim() const;
		// - SysSize returns the number of sites in the system.
		int SysSize() const;
		// Configuration manipulation functions:
		// - FermFZ::RandomCfg will generate a configuration with the specified start population.
		Config* RandomCfg() const;
		// - FermFZ::PropMove will move a fermion from one site to another, spin specific.
		Diff PropMove(Config* Cfg) const;
		// - Diff2Cfg inherited from GeneralFerm.
	};
}
#endif