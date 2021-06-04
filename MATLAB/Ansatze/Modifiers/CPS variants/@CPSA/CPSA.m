classdef CPSA < Modifier
    % CPSA - a Modifier subclass that modifies configuration amplitudes
    % using a Restricted Boltzmann Machine architecture of visible neurons
    % and hidden spins, parameterised as a correlator product state. A for
    % algebraic parameterisation.
    %   Original by X Fang, current version updated by M Pei.
    %   Modifier is the overarching class. CPS itself has subvariants with
    %   symmetries and projections built into them.
    
    % ---------------------------------
    % Format for CPSA Modifier object:
    % - CPSA.Nv = number of "visible" spins.
    % - CPSA.Nh = number of "hidden" spins.
    % - CPSA.Np = number of parameters in the ansatz = 2*Nv + 2*Nh + (2*Nv * 2*Nh).
    % - CPSA.a = (Nv x (VDim-1)) matrix - visible site vector elements.
    % - CPSA.b = (Nh x (HDim-1)) matrix - hidden site vector elements.
    % - CPSA.W = ((VDim-1) x (HDim-1) x Nv x Nh) array - hidden-visible coupling matrix elements.
    % - CPSA.HDim = 3 - hidden unit dimension.
    % - CPSA.VDim = 3 - visible unit dimension.
    % - CPSA.Ind0 = 1 - the fixed / zeroed element index for each correlator.
    % - CPSA.IndV = (VDim x 1) vector - translates v + Ind0 to a correlator index.
    % - CPSA.Theta = (Nh x (HDim-1)) matrix - effective angles.
    % - CPSA.VisInds = (Nv x 1) vector - a record of the current visible correlator indices.
    % ---------------------------------
    % Format for Update is a vector of new effective angles ThetaP and an
    % updated local index vector VisInds.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - (Nv*(VDim-1) x 1) for d/da.
    % Arranged ([i,1],...,[i,VDim],[i+1,1],...).
    % - (Nh*(HDim-1) x 1) for d/db.
    % Arranged ([i,1],...,[i,HDim],[i+1,1],...).
    % - ((Nh*Nv*(VDim-1)*(HDim-1)) x 1) for d/dW.
    % Arranged
    % ([v(i),h(j),i,j],[v(i)+1,h(j),i,j],...,[v(i),h(j)+1,i,j],
    %           ...,[v(i+1),h(j),i+1,j],...[v(i),h(j+1),i,j+1],...).
    % ---------------------------------
    
    properties % Default to one visible, one hidden plus state with no input.
        VFlag = 1; % Flag for whether to vary the parameters specified in this modifier. Can optimise.
    end
    
    properties (SetAccess = protected)
        Nv = 1; % Number of visible sites.
        Nh = 1; % Number of hidden sites.
        Np = 8; % Number of parameters. Np = Nv * 2 + Nh * 2 + Nv * Nh * 4
        a = 0; % Visible site correlator terms, Nv x 2 matrix.
        b = 0; % Hidden site correlator terms, Nh x 2 matrix.
        W = 0; % Hidden-visible correlator terms, Nv x Nh cell array of 2x2 matrices.
        HDim = 2; % Hidden site dimensionality - default to 2.
        VDim = 2; % Visible site dimensionality - determined by Hilbert.
        Graph % Details connectivity of lattice - used to include symmetries.
    end
    
    properties (Hidden)
        Theta = 0; % Effective hidden angles, Nh x 2 matrix.
        VisInds = 0; % A log of the configuration is included to aid ratio calculations.
        ParamCap = 100; % Parameter cap to mitigate effects of erroneous parameter changes.
        ParamMin = 1e-6; % Minimum magnitude permitted for parameters.
        OptInds = zeros(8,1); % Individual parameter flags for variational purposes.
    end
    
    properties (Hidden, SetAccess = protected)
        FullCfg = @FullSpinCfg; % Default case assumes spin Hilbert.
        Ind0 = 1; % Sets the 'zero' of the visible / hidden sites.
        IndV; % Translates n + Ind0 to correlator indices.
    end
    
    methods
        function obj = CPSA(Hilbert,Graph,Params,VFlag)
            % Graph necessary for CPS subvariants as second argument.
            if nargin < 4 % Assume variational if no VFlag specified. What is nargin
                obj.VFlag = 1;
            elseif nargin == 4
                obj.VFlag = VFlag;
            end
            if strcmp(Hilbert.Type,'Ferm')
                error('This CPS architecture is incompatible with fermionic systems.')
            else
                if Hilbert.d < 3
                    error('The CPS Modifier is intended for systems with visible dimension 3 or greater.');
                else
                    obj.Nv = Hilbert.N; obj.VDim = Hilbert.d;
                    if strcmp(Hilbert.Type,'Bose')
                        obj.FullCfg = @FullBoseCfg; obj.Ind0 = 1;
                        obj.IndV = 0:(obj.VDim-1);
                        % Set correlator elements to 1 for n = 0.
                    elseif strcmp(Hilbert.Type,'Spin')
                        obj.FullCfg = @FullSpinCfg; 
                        if mod(HilbertObj.d,2) % If d is odd, set middle value (0) to 1. 
                            obj.IndV = [1:((obj.VDim-1)/2) 0 ((obj.VDim+1)/2):(obj.VDim-1)];
                            obj.Ind0 = (obj.VDim+1)/2;
                        else % Set lowest spin value to 1.
                            obj.IndV = 0:(obj.VDim-1);
                            obj.Ind0 = (obj.VDim+1)/2;
                        end
                    end
                end
            end
            if isfield(Params,'Alpha')
                obj.Nh = Params.Alpha * Hilbert.N;
            else
                obj.Nh = Params.Nh;
            end
            obj.Graph = Graph;
            obj = RandomInitPsiCPS(obj,Params);
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateCPS(obj,dP);
        end
        
        % PsiCfgUpdate: Update Modifier configuration information inside
        % Update.
        function obj = PsiCfgUpdate(obj,Update)
            obj.Theta = Update.Theta;
            obj.VisInds = Update.VisInds;
        end
        
        % PrepPsi: Initialise Modifier configuration information given a
        % starting Cfg.
        function obj = PrepPsi(obj,Cfg)
            obj = PrepPsiCPSA(obj,Cfg);
        end
        
        % PsiGenerate: Generate full normalised NQS amplitudes for a given
        % set of basis states.
        function [Psi] = PsiGenerate(obj,Basis)
            Psi = PsiGenerateCPSA(obj,Basis);
        end
        
        % PsiRatio: Ratio of amplitudes for two configurations separated by
        % Diff.
        function [Ratio,Update] = PsiRatio(obj,Diff)
            % Set up new index vector IndsP.
            VisIndsP = obj.VisInds;
            for d = 1:Diff.num
                VisIndsP(Diff.pos(d)) = obj.VisInds(Diff.pos(d)) + Diff.val(d);
            end
            ThetaP = obj.Theta;
            for h = 1:obj.Nh
                for k = 1:Diff.num
                    % Apply changes to each Theta vector.
                    for hd = 1:(obj.HDim-1)
                        if obj.IndV(obj.VisInds(Diff.pos(k))) ~= 0
                            ThetaP(h,hd) = ThetaP(h,hd) / obj.W(obj.IndV(obj.VisInds(Diff.pos(k))),hd,Diff.pos(k),h);
                        end
                        if obj.IndV(VisIndsP(Diff.pos(k))) ~= 0
                            ThetaP(h,hd) = ThetaP(h,hd) * obj.W(obj.IndV(VisIndsP(Diff.pos(k))),hd,Diff.pos(k),h);
                        end
                    end
                end
            end
            Ratio = 1;
            % Apply visible correlator part of the ratio.
            for k = 1:Diff.num
                if obj.IndV(obj.VisInds(Diff.pos(k))) ~= 0
                    Ratio = Ratio / obj.a(Diff.pos(k),obj.IndV(obj.VisInds(Diff.pos(k))));
                end
                if obj.IndV(VisIndsP(Diff.pos(k))) ~= 0
                    Ratio = Ratio * obj.a(Diff.pos(k),obj.IndV(VisIndsP(Diff.pos(k))));
                end
            end
            % Perform the trace over the hidden dimension and take the product of trace
            % ratios.
            Ratio = Ratio * prod((1+sum(ThetaP,2))./(1+sum(obj.Theta,2)));
            if isnan(Ratio) || isinf(Ratio)
                Ratio = 0;
            end
            Update.Theta = ThetaP;
            Update.VisInds = VisIndsP;
        end
        
        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)
            Cfg_Inds = obj.FullCfg(Cfg) + obj.Ind0;
            dLogp = zeros(obj.Np,1);
            % Calculate all dTheta terms beforehand to save time.
            dTheta = ones(obj.Nh,1) ./ (1 + sum(obj.Theta,2));
            % Calculate d/da.
            for v = 1:obj.Nv
                if (obj.IndV(Cfg_Inds(v)) ~= 0)
                    PInd = (v-1)*(obj.VDim-1) + obj.IndV(Cfg_Inds(v));
                    if obj.OptInds(PInd)~=0
                        dLogp(PInd) = 1/obj.a(v,obj.IndV(Cfg_Inds(v)));
                    end
                    % Calculate d/dW.
                    for h = 1:obj.Nh
                        for hd = 1:(obj.HDim-1)
                            PInd = obj.IndV(Cfg_Inds(v)) + (obj.VDim-1)*(hd - 1 + obj.Nv + (obj.HDim-1)*...
                                (v - 1 + obj.Nv*(h-1))) + (obj.HDim-1)*obj.Nh;
                            if obj.OptInds(PInd) ~= 0
                                dLogp(PInd) = dTheta(h) / obj.W(obj.IndV(Cfg_Inds(v)),hd,v,h);
                            end
                        end
                    end
                end
            end
            % Calculate d/db.
            for h = 1:obj.Nh
                for hd = 1:(obj.HDim-1)
                    PInd = obj.Nv*(obj.VDim-1) + (h-1)*(obj.HDim-1) + hd;
                    if obj.OptInds(PInd)~=0
                        dLogp(PInd) = dTheta(h)/obj.b(h,hd);
                    end
                end
            end
            % Do some forward error prevention for NaN or Inf elements by zeroing them:
            dLogp(isnan(dLogp)) = 0;
            dLogp(isinf(dLogp)) = 0;
        end
        
        % ChangeInds: change the gauged out value of the visible units.
        function [obj] = ChangeInds(obj,NewInds)
            if length(NewInds) ~= obj.VDim
                error(['New vector must have ' num2str(obj.VDim) ' values.']);
            end
            NewInds = reshape(NewInds,numel(NewInds),1); MFlag = 0;
            for p = 0:(obj.VDim-1)
                if isempty(find(NewInds==p))
                    MFlag = 1;
                    disp(['Missing entry for p = ' num2str(p) '.']);
                end
            end
            if MFlag == 1
                error('Input index vector has missing entries.');
            end
            obj.IndV = NewInds;
        end
        
        % ParamList: outputs an Np x 1 vector of parameter values.
        function [Params] = ParamList(obj)
            Params = ParamListCPS(obj);
        end
        
        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
            obj = ParamLoadCPS(obj,P);
        end
        
        % PropertyList: Output a struct with the relevant properties as
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'CPSA';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds;
            Properties.Nv = obj.Nv; Properties.Nh = obj.Nh; Properties.ParamMin = obj.ParamMin;
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end
end
