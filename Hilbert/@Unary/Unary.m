classdef Unary < Hilbert
    % Unary - object containing details of a spin configuration Hilbert
    % space that performs unary encoding.
    %   Hilbert is overarching class.
    
    properties (SetAccess = protected)
        Type = 'Unary';
        N = 1; % Default to single site if no input.
        S = 1/2; % Default to spin-1/2 if no input.
        d = 2; % Dimension of variable to be unary encoded.
        Sector = []; % Empty sector - assume all Sz sectors permitted if no input.
    end
    
    properties (Hidden, SetAccess = protected)
        PropMoveFunc = @FlipUnaryCfg; % Default proposed configuration move.
        FullCfgFunc = @FullSpinCfg; % Default configuration vector for reference state.
        RandomCfgFunc = @RandomUnarySpin; % Default random configuration generating function.
        Diff2CfgFunc = @Diff2CfgSpinH; % Default difference to configuration conversion function.
        CParams = 1; % Default random configuration parameters.
    end
    
    methods
        % Constructor for subclass:
        function obj = Unary(N0,D)         
            % Number of unary sites is given by N0*D.            
            N = N0*D;
            CParams.N = N; CParams.d = D; CParams.N0 = N0;
            obj.N = N; obj.d = D; obj.CParams = CParams;
        end
        
        % RandomCfg: Generate Cfg struct subject to Hilbert parameters.
        function [Cfg] = RandomCfg(obj)
            [Cfg] = obj.RandomCfgFunc(obj.CParams);
        end
        
        % PropMove: Given a Cfg struct, generate a proposed CfgP and
        % corresponding difference struct Diff.
        function [Diff,CfgP] = PropMove(obj,Cfg)
            [Diff,CfgP] = obj.PropMoveFunc(Cfg);
        end
        
        % FullCfg: Given a Cfg struct, create a vector representation of
        % the configuration.
        function [Cfg_vec] = FullCfg(obj,Cfg)
            [Cfg_vec] = obj.FullCfgFunc(Cfg);
        end
        
        % Diff2Cfg: Given a Cfg and Diff struct, combine to create a new
        % CfgP struct.
        function [CfgP] = Diff2Cfg(obj,Diff,Cfg)
            [CfgP] = obj.Diff2CfgFunc(Diff,Cfg);
        end
        
        % PropertyList: Output a struct with the relevant properties as
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.N = obj.N; Properties.S = obj.S;
            Properties.Sector = obj.Sector;
        end
    end
    
end