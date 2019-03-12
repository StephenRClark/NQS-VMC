classdef HeisenbergSH < Hamiltonian
    % HeisenbergSH - an object describing a spin-half Heisenberg
    % Hamiltonian with nearest neighbour coupling. Also encompasses XXZ and
    % XYZ Hamiltonians, though not longer range exchange interactions.
    %   Hamiltonian is overarching class. Currently only compatible with
    %   Spin Hilbert spaces.
    
    properties
        HParams % List of parameters that characterise the Hamiltonian.
        Hilbert % Contains configuration information for the Markov chain.
        Operator % List of Operator objects that compose the Hamiltonian.
        TFuncs % Cell list of time-dependent functions that modify the HParams. Can be left empty.
        % Note that each Operator needs its own Graph that specifies their
        % connectivity on the lattice. 
        % For S = 1/2 Heisenberg, HParams contains the coupling constants
        % between spins. Sets energy scale. Maximum 3 entries.
    end
    
    properties (Hidden)
        HParams0 % Original input Hamiltonian parameters - used as a reference when time evolving.
    end
    
    methods
        % Constructor for the HeisenbergSH object:
        function obj = HeisenbergSH(Hilbert,Graph,J,GTFlag,TFuncs)
            % GTFlag specifies if bipartite gauge transform to change signs
            % of SxSx and SySy terms is desired.
            obj.Hilbert = Hilbert;
            if nargin == 2 % Assume isotropic, antiferromagnetic and with gauge transform if no J specified.
                obj.HParams = 1; GTFlag = 1;
            elseif nargin == 3 % Default to gauge transform active.
                obj.HParams = J; GTFlag = 1;
            elseif nargin == 4
                obj.HParams = J;
            end
            obj.HParams0 = obj.HParams;
            if nargin == 5
                if numel(TFuncs)~=numel(obj.HParams)
                    error('Mismatch between number of parameters and parameter evolution functions.')
                else
                    obj.TFuncs = TFuncs;
                end
            else
                obj.TFuncs = cell(numel(obj.HParams),1); % Empty cells.
            end
            % Assign Bonds for nearest neighbour lookup table.
            if numel(obj.HParams) == 1 % Isotropic case.
                % Easier to pull from a single function than a bunch of
                % separate function calls.
                if GTFlag == 1
                    obj.Operator{1} = Operator2S(Hilbert,Graph,@SiSj_GT_OpMatEls);
                else
                    obj.Operator{1} = Operator2S(Hilbert,Graph,@SiSj_OpMatEls);
                end
            elseif numel(obj.HParams) == 2 % Two inputs is taken to mean XXZ.
                % Rather than separate SxSx and SySy calls, takes matrix
                % elements from S+S- + S-S+.
                obj.Operator{1} = Operator2S(Hilbert,Graph,@SpmSmp_OpMatEls);
                obj.Operator{2} = OperatorDg(Hilbert,Graph,@SzSz_CfgVal);
                if GTFlag == 1 % J terms are listed in order (xy,z)
                    obj.HParams(1) = -obj.HParams(1);
                end
                    
            elseif numel(obj.HParams) == 3 % Three inputs is taken to mean XYZ.
                % Pulls matrix elements from individual SS operators.
                obj.Operator{1} = Operator2S(Hilbert,Graph,@SxSx_OpMatEls);
                obj.Operator{2} = Operator2S(Hilbert,Graph,@SySy_OpMatEls);
                obj.Operator{3} = OperatorDg(Hilbert,Graph,@SzSz_CfgVal);
                if GTFlag == 1 % J terms are listed in order (x,y,z)
                    obj.HParams(1) = -obj.HParams(1);
                    obj.HParams(2) = -obj.HParams(2);
                end
            end
            if strcmp(Hilbert.Type,'Spin') == 0
                error('Heisenberg model is currently only implemented for spin systems.')
            end            
        end
        
        % Modify HParams according to any TFuncs specified at
        % initialisation. Only called during time-dependent optimisations.
        function obj = TimeEvolveH(obj,time)
            for p = 1:numel(obj.HParams)
                if isempty(obj.TFuncs{p}) == 0
                    obj.HParams(p) = obj.TFuncs{p}(obj.HParams0(p),time);
                end
            end
        end
        
    end
    
    methods (Static)
        % Generate the matrix elements and the differences created by the
        % action of the Hamiltonian on a configuration Cfg.
        function [Diff,HMatEls] = HamMatEls(obj,Cfg)
            for j = 1:numel(obj.HParams)
                if j == 1
                    [Diff,OpMatEls] = obj.Operator{j}.CorrMatEls(obj.Hilbert,Cfg,obj.Operator{j}.Graph.Bonds);
                    HMatEls = OpMatEls * obj.HParams(j);
                else
                    [DiffT,OpMatElsT] = obj.Operator{j}.CorrMatEls(obj.Hilbert,Cfg,obj.Operator{j}.Graph.Bonds);
                    HMatEls = [HMatEls; OpMatElsT*obj.HParams(j)];
                    Diff = [Diff; DiffT];
                end
            end
        end
        
        % Sample the individual operators composing the Hamiltonian and
        % output the local energy EnLoc.
        function [EnLoc] = EnergySample(obj,Cfg,Ansatz)
            EnLoc = 0;
            for j = 1:numel(obj.HParams)
                % Call in turn each Operator's LocalSample function over
                % their associated Graph's Bonds, then weight by J
                % contribution.
                EnLoc = EnLoc + obj.HParams(j)*obj.Operator{j}.LocalSample(Cfg,0,0,Ansatz);
            end
        end
    end    
end

