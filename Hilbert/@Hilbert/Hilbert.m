classdef (Abstract) Hilbert
    % Hilbert - object containing details of Hilbert space in which the
    % desired wavefunction resides. Encompasses configuration functions.
    %   Hilbert is overarching class - subclasses will define specific
    %   functions for the methods listed below.
    
    properties (Abstract, SetAccess = protected)
        Type % Field to identify if configurations are bosonic, fermionic or spin.
        N % Number of sites.
        d % Single site Hilbert space dimension. In bosonic case, can be used as limiter on occupation.
        Sector % Defines permitted sector of Hilbert space - useful for fixed N or Sz calculations.
    end
    
    properties (Abstract, Hidden, SetAccess = protected)
        PropMoveFunc % Propose a configuration move subject to Hilbert space restrictions.
        FullCfgFunc % Create a vector output from a Cfg struct.
        RandomCfgFunc % Generate a random Cfg struct subject to Hilbert constraints.
        Diff2CfgFunc % Convert a Diff and Cfg into new CfgP.
        CParams % Parameters for generating a random Cfg struct.
    end
    
    methods (Abstract)
        % RandomCfg: Generate Cfg struct subject to Hilbert parameters.
        [Cfg] = RandomCfg(obj)
        
        % PropMove: Given a Cfg struct, generate a proposed CfgP and
        % corresponding difference struct Diff.
        [Diff,CfgP] = PropMove(obj,Cfg)
       
        % FullCfg: Given a Cfg struct, create a vector representation of
        % the configuration.
        [Cfg_vec] = FullCfg(obj,Cfg)
        
        % Diff2Cfg: Given a Cfg and Diff struct, combine to create a new
        % CfgP struct.
        [CfgP] = Diff2Cfg(obj,Diff,Cfg)
        
    end
end