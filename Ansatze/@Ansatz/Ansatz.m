classdef Ansatz
    % Ansatz - class of object containing variational wavefunction
    % parameters, reference states, projectors and methods
    %   Ansatz is overarching class - subclasses will define specific
    %   functions for the methods listed below.
    
    properties % (Abstract)
        Reference % Reference state - options listed below:
        % Plus (equal superposition) / BECR (Bose condensate) /
        % SDet (Slater determinant) / Pfaf (Pairing amplitude Pfaffian)
        Modifier % Amplitude modifier - options listed below:
        % None (no amplitude modification) / NQS (RBM architecture) /
        % Jast (Jastrow) / Gutz (Gutzwiller)
        Graph % Details connectivity of sites.
        Hilbert % Details Hilbert space in which Ansatz is represented.
        Var % Counter for number of variational components in Ansatz.
        NpTotal % Total number of variational parameters in Ansatz.
    end
    
    methods
        % General constructor for Ansatz class:
        
        function obj = Ansatz(Ref,Mod,Hilbert,Graph,Params)
            % Specific input scheme required:
            % Ref - 1 x 2 cell with name of reference state as 1st entry,
            % then flag with 0 for fixed or 1 for variational as 2nd.
            % Mod - N_mod x 2 cell with similar convention.
            % Hilbert - generate a Hilbert object which the Ansatz draws
            % configuration information from.
            % Graph - details connectivity of the lattice on which the
            % model is applied.
            % Params - struct containing necessary information to
            % initialise references and modifiers.
            
            % Initialise Reference object according to input.
            if strcmp(Ref{1},'Plus')
                obj.Reference = Plus;
                if Ref{2} == 1
                    disp(['Equal superposition of all states requested as reference state.'
                        ' Variational flag ignored as reference cannot be made variational.'])
                end
            elseif strcmp(Ref{1},'BECR')
                if Ref{2} == 1
                    disp('Variational bosonic references are not currently supported.')
                end
                if strcmp(Hilbert.Type,'Bose')
                    obj.Reference = BECR(Hilbert,Params);
                else
                    error('Bosonic reference requested for a non-bosonic Hilbert space.')
                end
            elseif strcmp(Ref{1},'SDet') || strcmp(Ref{1},'Pfaf')
                if strcmp(Hilbert.Type,'Bose')
                    error(['Fermionic reference requested for a bosonic Hilbert space. '
                        ' Advise using the "Bose" reference state instead.']);
                elseif strcmp(Hilbert.Type,'Spin')
                    disp(['Fermionic reference requested for a spin Hilbert space. '
                        'Ensure that "Hilbert" object is populated properly with N fermions.']);
                end
                if strcmp(Ref{1},'SDet')
                    obj.Reference = SDet(Hilbert,Params,Ref{2});
                elseif strcmp(Ref{1},'Pfaf')
                    % Pfaf can be initialised to respect symmetries in
                    % Graph, so will pass along to Pfaf.
                    obj.Reference = Pfaf(Hilbert,Graph,Params,Ref{2});
                end
            else
                error(['Requested reference does not match any available state. '
                    'Valid references are "Plus","Bose","SDet","Pfaf".'])
            end
            if Ref{2} == 1
                NpTotal = Reference.Np;
            else
                NpTotal = 0;
            end
            
            % Initialise Modifier objects in order.
            Modifier = cell(size(Mod,1),1); Var = Ref{2};
            for m = 1:size(Mod,1)
                % NQS Modifiers and subtypes.
                if strcmp(Mod{m,1},'NQS')
                    Modifier{m} = NQS(Hilbert,Graph,Params,Mod{m,2});
                elseif strcmp(Mod{m,1},'NQSNH')
                    % Number-like hidden units - only used for bosonic
                    % systems.
                    if strcmp(Hilbert.Type,'Bose')
                        Modifier{m} = NQSNH(Hilbert,Graph,Params,Mod{m,2});
                    else
                        Modifier{m} = NQS(Hilbert,Graph,Params,Mod{m,2});
                    end
                elseif strcmp(Mod{m,1},'NQSTI')
                    % Translation invariant NQS - translations specified in
                    % Graph.
                    Modifier{m} = NQSTI(Hilbert,Graph,Params,Mod{m,2});
                elseif strcmp(Mod{m,1},'NQSTIRI')
                    % Rotation invariance - only applied for 2D currently.
                    if numel(Graph.Dim) ~= 2 || strcmp(Graph.Type,'HypCub') == 0
                        error('Rotation invariance not currently supported for this Graph type.')
                    else
                        % Function that adds rotated versions of already
                        % present Bonds.
                        [Graph] = RotateBonds(Graph);
                        Modifier{m} = NQSTI(Hilbert,Graph,Params,Mod{m,2});
                    end
                elseif strcmp(Mod{m,1},'NQSNHTI')
                    % Translation invariance with spin symmetry and
                    % number-like hidden units.
                    if strcmp(Hilbert.Type,'Bose')
                        Modifier{m} = NQSNHTI(Hilbert,Graph,Params,Mod{m,2});
                    else
                        Modifier{m} = NQSTI(Hilbert,Graph,Params,Mod{m,2});
                    end
                elseif strcmp(Mod{m,1},'NQSTISS')
                    % Translation invariance with spin symmetry
                    if strcmp(Hilbert.Type,'Ferm') && (Hilbert.Sector(1) == Hilbert.Sector(2))
                        Modifier{m} = NQSTISS(Hilbert,Graph,Params,Mod{m,2});
                    else
                        Modifier{m} = NQSTI(Hilbert,Graph,Params,Mod{m,2});
                    end
                elseif strcmp(Mod{m,1},'NQSTIDDJ')
                    % Translation invariance with density-density like
                    % coupling terms.
                    if strcmp(Hilbert.Type,'Ferm') && (Hilbert.Sector(1) == Hilbert.Sector(2))
                        Modifier{m} = NQSTIDDJ(Hilbert,Graph,Params,Mod{m,2});
                    else
                        Modifier{m} = NQSTI(Hilbert,Graph,Params,Mod{m,2});
                    end
                    % Jastrow Modifiers and subtypes.
                elseif strcmp(Mod{m,1},'Jast')
                    Modifier{m} = Jast(Hilbert,Graph,Params,Mod{m,2});
                elseif strcmp(Mod{m,1},'JastTI')
                    % Translation invariant Jastrow factors - requires
                    % Graph.
                    Modifier{m} = JastTI(Hilbert,Graph,Params,Mod{m,2});
                elseif strcmp(Mod{m,1},'Gutz')
                    if strcmp(Hilbert.Type,'Ferm') == 0
                        error('Gutzwiller Modifier is only implemented for fermionic systems.')
                    else
                        Modifier{m} = Gutz(Hilbert,Params,Mod{m,2});
                    end
                elseif strcmp(Mod{m,1},'None')
                    Modifier{m} = None;
                end
                Var = Var + Mod{m,2};
                if Mod{m,2} == 1
                    NpTotal = NpTotal + Modifier{m}.Np;
                end
            end
            obj.Var = Var; obj.NpTotal = NpTotal; obj.Modifier = Modifier;            
            obj.Graph = Graph; obj.Hilbert = Hilbert;
        end
    end
    
    methods
        % Update Ansatz variational parameters according to changes dP.
        function [obj] = PsiUpdate(obj,dP)
            % Assume dP is a single column vector containing all the
            % parameters for each section.
            if obj.Reference.VFlag == 1
                P = obj.Reference.Np; % Number of entries in dP relevant to Reference.
                [obj.Reference] = obj.Reference.PsiUpdate( obj.Graph, dP(1:P) );
            else
                P = 0;
            end
            for m = 1:numel(obj.Modifier)
                if obj.Modifier{m}.VFlag == 1
                    p = obj.Modifier{m}.Np; % Number of entries in dP relative to this Modifier.
                    [obj.Modifier{m}] = obj.Modifier{m}.PsiUpdate( obj.Graph, dP((1:p)+P) );
                    P = P + p;
                end
            end
        end
        
        % Initialise Ansatz configuration values given a starting Cfg.
        function [obj] = PrepPsi(obj,Cfg)
            [obj.Reference] = PrepPsi(obj.Reference,obj.Hilbert,Cfg);
            for m = 1:numel(obj.Modifier)
                [obj.Modifier{m}] = PrepPsi(obj.Modifier{m},obj.Hilbert,Cfg);
            end
        end
        
        % Update Ansatz configuration values according to Update.
        function [obj] = PsiCfgUpdate(obj,Update)
            % Isolate update to reference state values.
            UpdateRf = Update{1}; [obj.Reference] = PsiCfgUpdate(obj.Reference,UpdateRf);
            % Incorporate contributions from projectors.
            for m = 1:numel(obj.Modifier)
                [obj.Modifier{m}] = PsiCfgUpdate(obj.Modifier{m},Update{m+1});
            end
        end
        
        % Activate / deactivate variational flag for Reference.
        function [obj] = VarFlagRef(obj,Flag)
            if Flag ~= 0 && Flag ~= 1
                error('Flag value should be 0 or 1 to deactivate / activate variational flag respectively.')
            end
            if strcmp(obj.Reference.Type,'Plus')
                error('Plus reference is locked - cannot change variational flag.')
            elseif strcmp(obj.Reference.Type,'Bose')
                error('Bose reference does not currently support variational versions - variational flag is locked.')
            elseif strcmp(obj.Reference.Type,'SDet') || strcmp(obj.Reference.Type,'Pfaf')
                % Editting variational flag.
                PFlag = obj.Reference.VFlag; obj.Reference.VFlag = Flag;
                % Updating number of variational components.
                Var0 = obj.Var - PFlag; obj.Var = Var0 + Flag;
                % Updating total number of variational parameters.
                NpTotal0 = obj.NpTotal - (PFlag * obj.Reference.Np);
                obj.NpTotal = NpTotal0 + (Flag * obj.Reference.Np);
                if Flag == 0
                    disp('Variational optimisation of reference has been disabled.')
                elseif Flag == 1
                    disp('Variational optimisation of reference has been enabled.')
                end
            end
        end
        
        % Activate / deactivate variational flag for Modifier specified by Mnum.
        function [obj] = VarFlagMod(obj,Mnum,Flag)
            if Flag ~= 0 && Flag ~= 1
                error('Flag value should be 0 or 1 to deactivate / activate variational flag respectively.')
            end
            if (Mnum > numel(obj.Modifier)) || mod(Mnum,1) ~= 0 || Mnum < 1
                error('Designated modifier does not exist.')
            end
            % Editting variational flag.
            PFlag = obj.Modifier{Mnum}.VFlag; obj.Modifier{Mnum}.VFlag = Flag;
            % Updating number of variational components.
            Var0 = obj.Var - PFlag; obj.Var = Var0 + Flag;
            % Updating total number of variational parameters.
            NpTotal0 = obj.NpTotal - (PFlag * obj.Modifier{Mnum}.Np);
            obj.NpTotal = NpTotal0 + (Flag * obj.Modifier{Mnum}.Np);
            if Flag == 0
                disp('Variational optimisation of modifier has been disabled.')
            elseif Flag == 1
                disp('Variational optimisation of modifier has been enabled.')
            end
            
        end
        
        % Exchange the existing Reference of an Ansatz object with a new
        % pre-prepared Reference.
        function obj = RefReplace(obj,NewRef)
            % Log the number of variational components and parameters when
            % old Reference is removed.
            Np0 = obj.NpTotal - (obj.Reference.Np * obj.Reference.VFlag);
            Var0 = obj.Var - obj.Reference.VFlag;
            % Check NewRef and Hilbert are compatible.
            HType = obj.Hilbert.Type; NRType = NewRef.Type;
            if strcmp(HType,'Ferm')
                if strcmp(NRType,'BECR')
                    error('New reference is only compatible with bosonic Hilbert spaces.')
                end
            elseif strcmp(HType,'Spin')
                if strcmp(NRType,'BECR')
                    error('New reference is only compatible with bosonic Hilbert spaces.')
                elseif strcmp(NRType,'Pfaf') || strcmp(NRType,'SDet')
                    error('New reference is only compatible with fermionic Hilbert spaces.')
                end
            elseif strcmp(HType,'Bose')
                if strcmp(NRType,'Pfaf') || strcmp(NRType,'SDet')
                    error('New reference is only compatible with fermionic Hilbert spaces.')
                end
            end
            obj.Reference = NewRef;
            obj.Var = Var0 + NewRef.VFlag; obj.NpTotal = Np0 + NewRef.VFlag * NewRef.Np;
        end
        
        % Exchange an existing Modifier of an Ansatz object with a new
        % pre-prepared Modifier.
        function obj = ModReplace(obj,NewMod,Nmod)
            % Nmod signifies which Modifier is to be replaced. If Nmod >
            % actual number of Modifiers, this will just append NewMod to
            % the list.
            if Nmod > numel(obj.Modifier)
                disp('Ansatz Modifier count is less than specified - appending new Modifier to list.')
                Np0 = obj.NpTotal; Var0 = obj.Var;
                obj.Modifier{numel(obj.Modifier)+1} = NewMod;
            else
                OldMod = obj.Modifier{Nmod};
                Np0 = obj.NpTotal - OldMod.Np; Var0 = obj.Var - OldMod.VFlag;
                obj.Modifier{Nmod} = NewMod;
            end
            obj.Var = Var0 + NewMod.VFlag; obj.NpTotal = Np0 + NewMod.Np;
        end
        
        % Remove an existing Modifier of an Ansatz object, with removed
        % Modifier as second output.
        function [obj,OldMod] = ModRemove(obj,Nmod)
            % Nmod signifies which Modifier is to be replaced. If Nmod >
            % actual number of Modifiers, this will just append NewMod to
            % the list.
            if Nmod > numel(obj.Modifier)
                error('Requested Modifier does not exist.')
            else
                OldMod = obj.Modifier{Nmod};
                obj.NpTotal = obj.NpTotal - OldMod.Np; obj.Var = obj.Var - OldMod.VFlag;
            end
            if numel(obj.Modifier) == 1
                % Replace with a None Modifier if last Modifier is removed.
                obj.Modifier{1} = None;
            else % Remove Modifier from list.
                obj.Modifier(Nmod) = [];
            end
        end
        
    end
    
    methods (Static)
        
        % Ratio between two configurations differing by Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            % Ratio from reference state:
            Update = cell(numel(obj.Modifier)+1,1);
            [RatioRf,Update{1}] = obj.Reference.PsiRatio(obj.Reference,Diff);
            RatioMd = zeros(numel(obj.Modifier),1);
            % Ratios from chosen Modifiers:
            for m = 1:numel(obj.Modifier)
                [RatioMd(m), Update{m+1}] = obj.Modifier{m}.PsiRatio(obj.Modifier{m},Diff);
            end
            % Combine ratios as one value and updates as a cell.
            Ratio = RatioRf * prod(RatioMd);
        end
        
        % Logarithmic derivative for the variational parameters in Ansatz.
        function [dLogp] = LogDeriv(obj,Cfg)
            dLogp = cell(obj.Var,1);
            if obj.Reference.VFlag == 1
                p = 1;
                dLogp{p} = obj.Reference.LogDeriv(obj.Reference,obj.Hilbert,obj.Graph,Cfg);
            else
                p = 0;
            end
            for m = 1:numel(obj.Modifier)
                if obj.Modifier{m}.VFlag == 1
                    p = p + 1;
                    [dLogpMod] = obj.Modifier{m}.LogDeriv(obj.Modifier{m},obj.Hilbert,obj.Graph,Cfg);
                    dLogp{p} = dLogpMod;
                end
            end
            dLogp = cell2mat(dLogp);
            % Ensure all subordinate LogDeriv functions output a Np x 1 vector.
        end
        
    end
    
end
