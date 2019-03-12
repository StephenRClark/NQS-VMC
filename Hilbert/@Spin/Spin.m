classdef Spin < Hilbert
    % Spin - object containing details of a spin configuration Hilbert
    % space.
    %   Hilbert is overarching class.
    
    properties
        Type = 'Spin';
        N = 1; % Default to single site if no input.
        S = 1/2; % Default to spin-1/2 if no input.
        d = 2; % Single site Hilbert space dimension.
        Sector = []; % Empty sector - assume all Sz sectors permitted if no input.
    end
    
    properties (Hidden)
        PropMove = @FlipSpinCfg; % Default proposed configuration move.
        FullCfgRef = @FullSpinCfg; % Default configuration vector for reference state.
        FullCfgMod = @FullSpinCfg; % Default configuration vector for amplitude modifier.
        RandomCfgFun = @RandomSpinZeroMag; % Default random configuration generating function.
        CParams = 1; % Default random configuration parameters.
    end
    
    methods
        % Constructor for subclass:
        function obj = Spin(N,S,Sector)
            % Currently only spin-1/2 and spin-1 are supported.
            if nargin == 1 % Default to spin-1/2 if no S defined.
                obj.S = 1/2;
                if mod(N,1)
                    error('Number of sites N should be integer.')
                else
                    obj.N = N; obj.Sector = [];
                end
            elseif nargin == 2 % Default to no Sz projection if no Sector input.
                obj.S = S; obj.N = N; obj.Sector = [];
            elseif nargin == 3
                obj.S = S; obj.N = N; obj.Sector = Sector;
            elseif nargin == 0
                obj.S = 1/2; obj.N = 1; obj.Sector = [];
            end
            if obj.S ~= 1/2 && obj.S ~= 1
                error('Only spin-1/2 and spin-1 are currently supported. Ensure the second argument is 1/2 or 1.')
            end
            obj.d = 2*obj.S + 1;            
            % Assign PropMove according to whether sector is defined or not.
            if isempty(obj.Sector)
                obj.PropMove = @FlipSpinCfg;
            else
                if obj.S == 1/2
                    obj.PropMove = @SwapSpinCfg;
                elseif obj.S == 1
                    obj.PropMove = @MoveSpin1Cfg;
                end
            end            
            % Assign RandomCfgFun and starting CParams according to S and
            % Sector respectively.
            if obj.S == 1/2
                obj.RandomCfgFun = @RandomSpinHFixedMag;
            elseif obj.S == 1
                obj.RandomCfgFun = @RandomSpin1FixedMag;
            end
            CParams.N = N;
            if isempty(obj.Sector)
                CParams.SzT = 0; % Default to starting at zero SzT if no Sector.
            else
                CParams.SzT = obj.Sector;
            end
            obj.CParams = CParams;
        end
    end
    
    methods (Static)
        % Call RandomCfg function according to sector.
        function [Cfg] = RandomCfg(obj)
            [Cfg] = obj.RandomCfgFun(obj.CParams);
        end
    end
end