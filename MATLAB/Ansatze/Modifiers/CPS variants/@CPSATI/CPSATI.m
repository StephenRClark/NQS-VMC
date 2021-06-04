classdef CPSATI < CPSA
    % CPS - a Modifier subclass that modifies configuration amplitudes
    % using a Restricted Boltzmann Machine architecture of visible neurons
    % and hidden spins, parameterised as a correlator product state. This
    % version features translation invariance.
    %   Original by X Fang, current version updated by M Pei.
    %   CPS is overarching class, which is itself a subclass of Modifier.
    
    % ---------------------------------
    % Format for CPSA Modifier object:
    % - CPSA.Nv = number of "visible" spins.
    % - CPSA.Nh = number of "hidden" spins.
    % - CPSA.Np = number of parameters in the ansatz = 2*Nv + 2*Nh + (2*Nv * 2*Nh).
    % - CPSA.a = (1 x (VDim-1)) vector - visible site vector elements.
    % - CPSA.b = (Alpha x (HDim-1)) matrix - hidden site vector elements.
    % - CPSA.W = ((VDim-1) x (HDim-1) x Nv x Alpha) array - hidden-visible coupling matrix elements.
    % - CPSA.HDim = 3 - this version features fixed hidden unit dimension.
    % - CPSA.VDim = 3 - this version is only compatible with Hilberts with dim = 3.
    % - CPSA.Ind0 = 1 - the fixed / zeroed element index for each correlator.
    % - CPSA.IndV = (VDim x 1) vector - translates v + Ind0 to a correlator index.
    % - CPSA.Theta = (Nh x (HDim-1)) matrix - effective angles.
    % - CPSA.VisInds = (Nv x 1) vector - a record of the current visible correlator indices.
    % Properties added with translation invariance:
    % - CPSA.Alpha = number of distinct hidden unit sets.
    % - CPSA.ati = (1 x (VDim-1)) vector - reduced parameter set for TI.
    % - CPSA.bti = (Alpha x (HDim-1)) matrix - reduced parameter set for TI.
    % - CPSA.Wti = ((VDim-1) x (HDim-1) x Nv x Alpha) array - reduced parameter set for TI.
    % ---------------------------------
    % Format for Update is a vector of new effective angles ThetaP and an
    % updated local index vector VisInds.
    % ---------------------------------
    % Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
    % - ((VDim-1) x 1) for d/da.
    % Arranged ([i,1],...,[i,VDim],[i+1,1],...).
    % - (Alpha*(HDim-1) x 1) for d/db.
    % Arranged ([i,1],...,[i,HDim],[i+1,1],...).
    % - ((Alpha*Nv*(VDim-1)*(HDim-1)) x 1) for d/dW.
    % Arranged
    % ([v(i),h(a),i,a],[v(i)+1,h(a),i,a],...,[v(i),h(a)+1,i,a],
    %           ...,[v(i+1),h(a),i+1,a],...[v(i),h(a+1),i,a+1],...).
    % ---------------------------------
    
    properties (SetAccess = protected)
        ati = 0; % Reduced TI parameter set for a, 1 x 2 vector.
        bti = 0; % Reduced TI parameter set for b, Alpha x 2 matrix.
        Wti = 0; % Reduced TI parameter set for W, Alpha x Nv cell array of 2x2 matrices.
        Alpha = 1;
    end
    
    methods
        function obj = CPSATI(Hilbert,Graph,Params,VFlag)
            obj@CPSA(Hilbert,Graph,Params,VFlag);
            obj = RandomInitPsiCPSTI(obj,Params);
        end
        
        % PsiUpdate: Update Modifier variational parameters according to
        % changes dP.
        function obj = PsiUpdate(obj,dP)
            obj = PsiUpdateCPSTI(obj,dP);
        end
        
        % PsiCfgUpdate inherited from CPSA.
        
        % PrepPsi inherited from CPSA.
        
        % PsiGenerate inherited from CPSA.
        
        % PsiRatio inherited from CPSA.
        
        % LogDeriv: Logarithmic derivative for the variational parameters
        % in Modifier.
        function [dLogp] = LogDeriv(obj,Cfg)
            GraphObj = obj.Graph; BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
            % translates by some combination of Graph.Lvecs.
            Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.
            Cfg_Inds = obj.FullCfg(Cfg) + obj.Ind0;
            dLogp = zeros(obj.Np,1);
            % Calculate d/da.
            for vd = 1:(obj.VDim-1)
                if obj.OptInds(vd) ~= 0
                    dLogp(vd) = sum(obj.IndV(Cfg_Inds)==vd)/obj.ati(vd);
                end
            end
            % Calculate all dTheta terms beforehand to save time.
            dTheta = obj.Theta ./ (1 + sum(obj.Theta,2));
            % Calculate d/db.
            for al = 1:obj.Alpha
                for hd = 1:(obj.HDim-1)
                    PInd = (obj.VDim-1) + (al-1)*(obj.HDim-1) + hd;
                    if obj.OptInds(PInd)~=0
                        dLogp(PInd) = sum(dTheta((1:Ntr)+(al-1)*Ntr))/obj.bti(al);
                    end
                    % Calculate d/dW.
                    for v = 1:obj.Nv
                        for bn = 1:Ntr
                            VInd = BondMap{bn}(v); HInd = bn + (al-1)*Ntr;
                            if (obj.IndV(Cfg_Inds(VInd)) ~= 0)
                                PInd = obj.IndV(Cfg_Inds(VInd)) + (obj.VDim-1)*(hd + (obj.HDim-1)*(...
                                    v - 1 + obj.Nv*(al-1))) + (obj.HDim-1)*obj.Alpha;
                                if obj.OptInds(PInd) ~= 0
                                    dLogp(PInd) = dLogp(PInd) + dTheta(HInd)/obj.Wti(obj.IndV(Cfg_Inds(VInd)),hd,VInd,al);
                                end
                            end
                        end
                    end
                end
            end
            % Do some forward error prevention for NaN or Inf elements by zeroing them:
            dLogp(isnan(dLogp)) = 0;
            dLogp(isinf(dLogp)) = 0;
        end
        
        % ParamList: outputs an Np x 1 vector of parameter values.
        function [Params] = ParamList(obj)
            Params = ParamListCPSTI(obj);
        end
        
        % ParamLoad: replaces parameters with the provided ones in vector P.
        function [obj] = ParamLoad(obj,P)
           obj = ParamLoadCPSTI(obj,P);
        end
        
        % PropertyList: Output a struct with the relevant properties as 
        % separate fields. Used for interfacing with C++ code.
        function [Properties] = PropertyList(obj)
            Properties.Type = 'CPSATI';
            Properties.Graph = obj.Graph.PropertyList; Properties.OptInds = obj.OptInds;
            Properties.Nv = obj.Nv; Properties.Alpha = obj.Alpha; 
            Properties.Params = obj.ParamList; Properties.ParamCap = obj.ParamCap;
        end
    end
end
