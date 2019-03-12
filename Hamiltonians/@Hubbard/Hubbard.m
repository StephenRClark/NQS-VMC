classdef Hubbard < Hamiltonian
    % Hubbard - an object describing a Hubbard Hamiltonian with nearest
    % neighbour hopping. Also includes options for longer range
    % density-density correlations.
    %   Hamiltonian is overarching class. Not compatible with Spin Hilbert
    %   spaces.
    
    properties
        HParams % List of parameters that characterise the Hamiltonian.
        TFuncs % Cell list of time-dependent functions that modify the HParams. Can be left empty.
        Hilbert % Contains configuration information for the Markov chain.
        Operator % List of Operator objects that compose the Hamiltonian.
        % Note that each Operator needs its own Graph that specifies their
        % connectivity on the lattice.
        % For Hubbard model, 1st entry is hopping amplitude(s) t, second is
        % on-site density-density cost U, all further entries are
        % density-density correlations with separations defined by Graphs.
    end
    
    properties (Hidden)
        HParams0 % Original input Hamiltonian parameters - used as a reference when time evolving.
    end
    
    methods
        % Constructor for the Hubbard object:
        function obj = Hubbard(Hilbert,Graphs,t,U,TFuncs)
            % U should be entered as a row vector.
            obj.Hilbert = Hilbert;
            obj.HParams = [t, U];
            if nargin == 5
                if numel(TFuncs)~=numel(obj.HParams)
                    error('Mismatch between number of parameters and parameter evolution functions.')
                else
                    obj.TFuncs = TFuncs;
                end
            else
                obj.TFuncs = cell(numel(obj.HParams),1); % Empty cells.
            end
            obj.HParams0 = obj.HParams;
            % Check number of Graphs inputted matches with number of
            % Operators needed.
            if numel(Graphs) < numel(t)+numel(U)
                error('Mismatch between number of Graphs and Operators.')
            end
            % Assign hopping operators to OpMatEls slots. Each hopping term
            % requires its own Graph to be defined.
            if strcmp(Hilbert.Type,'Ferm')
                for h = 1:numel(t)
                    obj.Operator{h} = Operator2S(Hilbert,Graphs(h),@CpCm_UD_OpMatEls);
                end
            elseif strcmp(Hilbert.Type,'Bose')
                for h = 1:numel(t)
                    obj.Operator{h} = Operator2S(Hilbert,Graphs(h),@BpBm_OpMatEls);
                end
            end
            % 
            if strcmp(Hilbert.Type,'Ferm')
                for h = 1:numel(U)
                    if h == 1 % First density-density term in Fermi case only counts double occupancies.
                        DDFunc = @DbOc_Ferm_CfgVal;
                    else
                        DDFunc = @NiNj_Ferm_CfgVal;
                    end
                    obj.Operator{h+numel(t)} = OperatorDg(Hilbert,Graphs(h+numel(t)),DDFunc);
                end
            elseif strcmp(Hilbert.Type,'Bose')
                for h = 1:numel(U)
                    obj.Operator{h+numel(t)} = OperatorDg(Hilbert,Graphs(h+numel(t)),@NiNj_Bose_CfgVal);
                end
            end
            if strcmp(Hilbert.Type,'Spin')
                error('Hubbard Hamiltonian is incompatible with spin systems.')
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
            % Calculate the matrix elements for the hopping terms.
            [Diff,CorrMatEls] = obj.Operator{1}.CorrMatEls(obj.Hilbert,Cfg,obj.Operator{1}.Graph);
            HMatEls = CorrMatEls * obj.HParams(1);
            % Calculate matrix elements for the density-density operator(s).
            for h = 2:numel(obj.HParams)
                OpInd = min(numel(obj.Operator),h);
                [DiffH,CorrMatElsH] = obj.Operator{OpInd}.CorrMatEls(obj.Hilbert,Cfg,obj.Operator{OpInd}.Graph);
                Diff = [Diff, DiffH]; HMatEls = [HMatEls; CorrMatElsH*obj.HParams(h)];
            end
        end
        
        % Sample the individual operators composing the Hamiltonian and
        % output the local energy EnLoc.
        function [EnLoc] = EnergySample(obj,Cfg,Ansatz)
            EnLoc = 0;
            for h = 1:numel(obj.HParams)
                % Call in turn each Operator's LocalSample function over
                % their associated Graph's Bonds, then weight by J
                % contribution.
                EnLoc = EnLoc + obj.HParams(h)*obj.Operator{h}.LocalSample(Cfg,0,0,Ansatz);
            end
        end
    end
    
end

