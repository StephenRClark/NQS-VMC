% --- General CPS logarithmic derivative function ---

% - Original by X Fang, updated by M Pei.

function dLogp = LogDerivCPSX(CPSObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the CPS ansatz, for a  configuration specifed
% by the structure Cfg.
% X indicates exponential parameterisation.
% ---------------------------------
% Format for CPSX Modifier object:
% - CPSX.Nv = number of "visible" spins.
% - CPSX.Nh = number of "hidden" spins.
% - CPSX.Np = number of parameters in the ansatz = 2*Nv + 2*Nh + (2*Nv * 2*Nh).
% - CPSX.a = (Nv x (VDim-1)) matrix - visible site vector elements.
% - CPSX.b = (Nh x (HDim-1)) matrix - hidden site vector elements.
% - CPSX.W = ((VDim-1) x (HDim-1) x Nv x Nh) array - hidden-visible coupling matrix elements.
% - CPSX.HDim = 3 - this version features fixed hidden unit dimension.
% - CPSX.VDim = 3 - this version is only compatible with Hilberts with dim = 3.
% - CPSX.Ind0 = 1 - the fixed / zeroed element index for each correlator.
% - CPSX.IndV = (VDim x 1) vector - translates v + Ind0 to a correlator index.
% - CPSX.Theta = (Nh x (HDim-1)) matrix - effective angles.
% - CPSX.VisInds = (Nv x 1) vector - a record of the current visible correlator indices.
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

Nv = CPSObj.Nv; VDim = CPSObj.VDim;
Nh = CPSObj.Nh; HDim = CPSObj.HDim;
IndV = CPSObj.IndV; Ind0 = CPSObj.Ind0;

Cfg_Inds = CPSObj.FullCfg(Cfg) + CPSObj.Ind0;
dLogp = zeros(CPSObj.Np,1);
% Calculate all dTheta terms beforehand to save time.
dTheta = exp(obj.Theta) ./ (1 + sum(exp(obj.Theta),2));
% Calculate d/da.
for v = 1:Nv
    if (CPSObj.IndV(Cfg_Inds(v)) ~= 0)
        PInd = (v-1)*(VDim-1) + IndV(Cfg_Inds(v));
        if CPSObj.OptInds(PInd)~=0
            dLogp(PInd) = 1;
        end
        % Calculate d/dW.
        for h = 1:Nh
            for hd = 1:(HDim-1)
                PInd = IndV(Cfg_Inds(v)) + (VDim-1)*(hd - 1 + Nv + (HDim-1)*...
                    (v - 1 + Nv*(h-1))) + (HDim-1)*Nh;
                if CPSObj.OptInds(PInd) ~= 0
                    dLogp(PInd) = dTheta(h,hd);
                end
            end
        end
    end
end
% Calculate d/db.
for h = 1:Nh
    for hd = 1:(HDim-1)
        PInd = Nv*(VDim-1) + (h-1)*(HDim-1) + hd;
        if CPSObj.OptInds(PInd)~=0
            dLogp(PInd) = dTheta(h,hd);
        end
    end
end
% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;
end
