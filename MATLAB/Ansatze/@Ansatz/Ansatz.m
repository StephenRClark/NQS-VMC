classdef Ansatz
    % Ansatz - class of object containing variational wavefunction
    % parameters, reference states, projectors and methods
    %   Ansatz is overarching class - subclasses will define specific
    %   functions for the methods listed below.
    
    properties (SetAccess = protected)
        Reference % Reference state - options listed below:
        % Plus (equal superposition) / BECR (Bose condensate) /
        % SDet (Slater determinant) / Pfaf (Pairing amplitude Pfaffian)
        Modifier % Amplitude modifier - options listed below:
        % None (no amplitude modification) / NQS (RBM architecture) /
        % Jast (Jastrow) / Gutz (Gutzwiller)
        Hilbert % Details Hilbert space in which Ansatz is represented.
        Var % Counter for number of variational components in Ansatz.
        NpTotal % Total number of variational parameters in Ansatz.
    end
    
    methods
        % General constructor for Ansatz class:        
        function obj = Ansatz(Ref,Mods,Hilbert)
            % To be cleaned up and replaced later - next version will
            % initialise with pre-made Reference and Modifiers.
            obj.Hilbert = Hilbert; % Common Hilbert space for Ref and Mod to use.
            if numel(Ref) > 1
                error('Only one Reference may be used per Ansatz.')
            end
            HType = Hilbert.Type; RType = Ref.Type; % Use to gauge compatibility.
            if strcmp(RType,'Plus') == 0
                if strcmp(RType,HType) == 0
                    error('Chosen Reference and Hilbert are incompatible.')
                end
            end
            obj.Reference = Ref; obj.Modifier = Mods;
            % Modifiers should be initialised with same Hilbert as input to
            % avoid issues arising - responsibility for this falls on the
            % user, though.
            Var = Ref.VFlag; NpTotal = Ref.Np * Ref.VFlag;
            for m = 1:numel(Mods)
                Var = Var + Mods{m}.VFlag;
                NpTotal = NpTotal + (Mods{m}.VFlag*Mods{m}.Np);
            end
            obj.Var = Var; obj.NpTotal = NpTotal;
        end
    end
    
    methods
        % PsiUpdate: Update Ansatz variational parameters according to
        % changes dP. Wrapper function that incorporates Reference and
        % Modifier versions of this function.
        function [obj] = PsiUpdate(obj,dP)
            % Assume dP is a single column vector containing all the
            % parameters for each section.
            if obj.Reference.VFlag == 1
                P = obj.Reference.Np; % Number of entries in dP relevant to Reference.
                [obj.Reference] = obj.Reference.PsiUpdate(dP(1:P));
            else
                P = 0;
            end
            for m = 1:numel(obj.Modifier)
                if obj.Modifier{m}.VFlag == 1
                    p = obj.Modifier{m}.Np; % Number of entries in dP relative to this Modifier.
                    [obj.Modifier{m}] = obj.Modifier{m}.PsiUpdate(dP((1:p)+P));
                    P = P + p;
                end
            end
        end
        
        % PsiUpdate: Update Ansatz variational parameters according to
        % changes dP. Wrapper function that incorporates Reference and
        % Modifier versions of this function.
        function [obj] = ParamLoad(obj,dP)
            % Assume dP is a single column vector containing all the
            % parameters for each section.
            if obj.Reference.VFlag == 1
                P = obj.Reference.Np; % Number of entries in dP relevant to Reference.
                [obj.Reference] = obj.Reference.ParamLoad(dP(1:P));
            else
                P = 0;
            end
            for m = 1:numel(obj.Modifier)
                if obj.Modifier{m}.VFlag == 1
                    p = obj.Modifier{m}.Np; % Number of entries in dP relative to this Modifier.
                    [obj.Modifier{m}] = obj.Modifier{m}.ParamLoad(dP((1:p)+P));
                    P = P + p;
                end
            end
        end
        
        % PrepPsi: Initialise Ansatz configuration values given a starting
        % Cfg. Wrapper function that incorporates Reference and Modifier
        % versions of this function.
        function [obj] = PrepPsi(obj,Cfg)
            [obj.Reference] = PrepPsi(obj.Reference,Cfg);
            for m = 1:numel(obj.Modifier)
                [obj.Modifier{m}] = PrepPsi(obj.Modifier{m},Cfg);
            end
        end
        
        % PsiCfgUpdate: update Ansatz configuration values according to
        % Update. Wrapper function that incorporates Reference and Modifier
        % versions of this function.
        function [obj] = PsiCfgUpdate(obj,Update)
            % Isolate update to reference state values.
            UpdateRf = Update{1}; [obj.Reference] = PsiCfgUpdate(obj.Reference,UpdateRf);
            % Incorporate contributions from projectors.
            for m = 1:numel(obj.Modifier)
                [obj.Modifier{m}] = PsiCfgUpdate(obj.Modifier{m},Update{m+1});
            end
        end
        
        % VarFlagRef: Activate / deactivate variational flag for Reference.
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
        
        % VarFlagMod: Activate / deactivate variational flag for Modifier
        % specified by Mnum.
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
                disp(['Variational optimisation of Modifier ' num2str(Mnum) ' has been disabled.'])
            elseif Flag == 1
                disp(['Variational optimisation of Modifier ' num2str(Mnum) ' has been enabled.'])
            end
            
        end
        
        % RefReplace: Exchange the existing Reference of an Ansatz object
        % with a new pre-prepared Reference.
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
        
        % ModReplace: Exchange an existing Modifier of an Ansatz object
        % with a new pre-prepared Modifier.
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
        
        % ModRemove: Remove an existing Modifier of an Ansatz object, with
        % removed Modifier as second output.
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
        
        % RndBatchSelect: alter OptInds of each component for random batch
        % optimisation.
        function [obj] = RndBatchSelect(obj,OFrac)
            if obj.Reference.VFlag == 1 % Only alter if Reference is being flagged for optimisation.
                obj.Reference = obj.Reference.RndBatchSelect(OFrac);
            end
            for m = 1:numel(obj.Modifier)
                if obj.Modifier{m}.VFlag == 1
                    obj.Modifier{m} = obj.Modifier{m}.RndBatchSelect(OFrac);
                end
            end
        end
        
        % ParamList; outputs an NpTotal x 1 vector of parameter values from
        % all actively optimisable constituents.
        function [Params] = ParamList(obj)
            Params = cell(obj.Var,1);
            p = 0;
            if obj.Reference.VFlag == 1
                p = p+1;
                Params{p} = obj.Reference.ParamList;
            end
            for m = 1:numel(obj.Modifier)
                if obj.Modifier{m}.VFlag == 1
                    p = p + 1;
                    Params{p} = obj.Modifier{m}.ParamList;
                end
            end
            Params = cell2mat(Params); 
            % Ensure all other ParamList functions list as Np x 1 or this will go horribly wrong.
        end
        
        % HilbertReplace: Replace Hilbert with new one. New Hilbert must
        % have the same Type, number of sites N and dimension d to avoid
        % compatibility issues.
        function [obj] = HilbertReplace(obj,NewHilbert)
            OldHilbert = obj.Hilbert;
            if (OldHilbert.N ~= NewHilbert.N)
                error('Hilbert spaces have differing number of sites N.');
            elseif strcmp(OldHilbert.Type,NewHilbert.Type) ==0
                error('Hilbert spaces describe different system types.');
            elseif (OldHilbert.d ~= NewHilbert.d)
                error('Hilbert spaces have different on-site dimension.');
            else
                obj.Hilbert = NewHilbert;
            end
        end
 
        % PsiRatio: Ratio between two configurations differing by Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            % Ratio from reference state:
            Update = cell(numel(obj.Modifier)+1,1);
            [Ratio,Update{1}] = obj.Reference.PsiRatio(Diff);
            % Ratios from chosen Modifiers:
            for m = 1:numel(obj.Modifier)
                [RatioMd, Update{m+1}] = obj.Modifier{m}.PsiRatio(Diff);
                Ratio = Ratio * RatioMd;
            end
        end
        
        % LogDeriv: Logarithmic derivative for the variational parameters in Ansatz.
        function [dLogp] = LogDeriv(obj,Cfg)
            dLogp = cell(obj.Var,1);
            if obj.Reference.VFlag == 1
                p = 1;
                dLogp{p} = obj.Reference.LogDeriv(Cfg);
            else
                p = 0;
            end
            for m = 1:numel(obj.Modifier)
                if obj.Modifier{m}.VFlag == 1
                    p = p + 1;
                    [dLogpMod] = obj.Modifier{m}.LogDeriv(Cfg);
                    dLogp{p} = dLogpMod;
                end
            end
            dLogp = cell2mat(dLogp);
            % Ensure all subordinate LogDeriv functions output a Np x 1 vector.
        end 
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Hilbert = obj.Hilbert.PropertyList;
            Properties.Reference = obj.Reference.PropertyList;
            Properties.Modifier = cell(numel(obj.Modifier),1);
            for o = 1:numel(obj.Modifier)
                Properties.Modifier{o} = obj.Modifier{o}.PropertyList;
            end
        end
    end
    
end
