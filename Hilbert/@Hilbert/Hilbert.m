classdef (Abstract) Hilbert
    % Hilbert - object containing details of Hilbert space in which the
    % desired wavefunction resides. Encompasses configuration functions.
    %   Hilbert is overarching class - subclasses will define specific
    %   functions for the methods listed below.
    
    properties (Abstract)
        Type % Field to identify if configurations are bosonic, fermionic or spin.
        N % Number of sites.
        d % Single site Hilbert space dimension. In bosonic case, can be used as limiter on occupation.
        Sector % Defines permitted sector of Hilbert space - useful for fixed N or Sz calculations.
    end
    
    properties (Abstract, Hidden)
        PropMove % Propose a configuration move subject to Hilbert space restrictions.
        FullCfgRef % Create a vector output from a Cfg struct for use with the reference.
        FullCfgMod % Create a vector output from a Cfg struct for use with the modifier.
        RandomCfgFun % Generate a random Cfg struct subject to Hilbert constraints.
        CParams % Parameters for generating a random Cfg struct.
    end
    
    methods (Static)
        [Cfg] = RandomCfg(obj)
    end
end

