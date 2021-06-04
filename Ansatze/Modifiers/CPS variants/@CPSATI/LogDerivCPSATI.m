% --- General CPS logarithmic derivative function ---

% - Original by X Fang, updated by M Pei.

function dLogp = LogDerivCPSATI(CPSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the CPS ansatz, for a  configuration specifed
% by the structure Cfg.
% A indicates algebraic parameterisation.
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

Nv = CPSObj.Nv; VDim = CPSObj.VDim; HDim = CPSObj.HDim;
IndV = CPSObj.IndV; Alpha = CPSObj.Alpha;
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; % Bond map detailing all possible distinct
% translates by some combination of Graph.Lvecs.
Ntr = numel(BondMap); % Multiplicity of each hidden unit layer depends on number of distinct translates.

Cfg_Inds = CPSObj.FullCfg(Cfg) + CPSObj.Ind0;
dLogp = zeros(CPSObj.Np,1);
% Calculate d/da.
for vd = 1:(VDim-1)
    if CPSObj.OptInds(vd) ~= 0
        dLogp(vd) = sum(IndV(Cfg_Inds)==vd)/CPSObj.ati(vd);
    end
end
% Calculate all dTheta terms beforehand to save time.
dTheta = CPSObj.Theta ./ (1 + sum(CPSObj.Theta,2));
% Calculate d/db.
for a = 1:Alpha
    for hd = 1:(HDim-1)
        PInd = (VDim-1) + (a-1)*(HDim-1) + hd;
        if CPSObj.OptInds(PInd)~=0
            dLogp(PInd) = sum(dTheta((1:Ntr)+(a-1)*Ntr))/CPSObj.bti(a);
        end
        % Calculate d/dW.
        for v = 1:Nv
            for b = 1:Ntr
                VInd = BondMap{b}(v); HInd = b + (a-1)*Ntr;
                if (CPSObj.IndV(Cfg_Inds(VInd)) ~= 0)
                    PInd = IndV(Cfg_Inds(VInd)) + (VDim-1)*(hd + (HDim-1)*...
                        (v - 1 + Nv*(a-1))) + (HDim-1)*Alpha;
                    if CPSObj.OptInds(PInd) ~= 0
                        dLogp(PInd) = dLogp(PInd) + dTheta(HInd,hd)/CPSObj.Wti(IndV(Cfg_Inds(VInd)),hd,VInd,a);
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