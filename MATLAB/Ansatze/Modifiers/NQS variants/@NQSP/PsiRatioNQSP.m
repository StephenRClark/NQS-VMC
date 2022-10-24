% --- General NQS wave function amplitude ratio function ---

function [Ratio,Update] = PsiRatioNQSP(NQSObj,Diff)
% This function computes the ratio Psi(CfgP)/Psi(Cfg) of amplitudes for
% a proposed spin-1/2 configuration CfgP and the current on Cfg, whose
% difference is stored in Diff.
% ---------------------------------
% Format for NQSP Modifier:
% - NQSP.Nv = number of "visible" units.
% - NQSP.Nh = number of "hidden" units.
% - NQSP.Np = number of parameters in the ansatz = (Nsl x VOrder) + (Alpha x
% HOrder) + (Nv x VOrder)(Alpha x HOrder)
% - NQSP.VDim = dimension of the visible units.
% - NQSP.HDim = dimension of the hidden units.
% - NQSP.VOrder = highest power of visible unit interactions. Max value VDim-1.
% - NQSP.HOrder = highest power of hidden unit interactions. Max value HDim-1.
% - NQSP.a = (Nv x VOrder) matrix - visible site biases.
% - NQSP.av = (Nsl x VOrder) matrix - visible bias parameters
% - NQSP.b = (Nh x HOrder) matrix - hidden site bias.
% - NQSP.bv = (Alpha x HOrder) matrix - hidden bias parameters.
% - NQSP.W = (Nh x Nv x HOrder x VOrder) array - hidden-visible coupling terms.
% - NQSP.Wm = (Alpha x Nv x HOrder x VOrder) array - hidden-visible coupling parameters
% - NQSP.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSP.Theta = (Nh x HOrder) matrix - effective angles by hidden order.
% - NQSP.VisVec = (Nv x 1) vector - visible occupancies.
% - NQSP.Rescale = flag for visible unit rescaling to [0 1] interval.
% ---------------------------------
% Format for Update is a struct with two fields:
% Update.Theta - matrix of new effective angles ThetaP.
% Update.VisVec - vector of new visible occupancies.
% ---------------------------------

% Make local copies to reduce notation.
HOrder = NQSObj.HOrder; HDim = NQSObj.HDim;
VOrder = NQSObj.VOrder; VDim = NQSObj.VDim; 
Rescale = NQSObj.Rescale;

VisVecP = NQSObj.VisVec(:); 
VisPow_shift = zeros(Diff.num,VOrder); % Shift in visible occupancy powers for each altered site.
VisPow_shift(:,1) = Diff.val(:)*((VDim-1)^(-Rescale));
% Calculate power shifts and visible bias contributions.
Ratio = exp(sum(a(Diff.pos(:),1).*VisPow_shift(:,1)));
for vo = 2:VOrder
    VisPow_shift(:,vo) = (VisVecP(Diff.pos(:)+VisPow_shift(:,1))).^vo - VisVecP(Diff.pos(:)).^vo;
    Ratio = Ratio * exp(sum(a(Diff.pos(:),vo).*VisPow_shift(:,vo)));
end
Theta_shift = zeros(NQSObj.Nh,HOrder); % Initialise effective angle shift.
% Only loop over the sites where there are differences:
for d=1:Diff.num
    for ho = 1:HOrder
        Theta_shift(:,ho) = Theta_shift(:,ho) + sum(reshape(...
            NQSObj.W(:,Diff.pos(d),ho,:),NQSObj.Nh,VOrder).*VisPow_shift(d,:),2);
    end
end
ThetaP = NQSObj.Theta + Theta_shift; % Update the effective angles for the proposed configuration.
% Calculate hidden unit contributions to ratio.
HidVals = zeros(HDim,HOrder);
HidVals(:,1) = (0:(HDim-1)).'*((HDim-1)^(-Rescale));
for ho = 2:HOrder
    HidVals(:,ho) = HidVals(:,1).^ho;
end
HidAmp = zeros(NQSObj.Nh,1); HidAmpP = zeros(NQSObj.Nh,1);
for h = 1:HDim
    HidAmp = HidAmp + exp(NQSObj.Theta*(HidVals(h,:).'));
    HidAmpP = HidAmpP + exp(ThetaP*(HidVals(h,:).'));
end
Ratio = Ratio * prod(HidAmpP./HidAmp);
VisVecP(Diff.pos(:)) = VisVecP(Diff.pos(:))+Diff.val(:)*((VDim-1)^(-Rescale));
% Add new VisVec and Theta to Update.
Update.Theta = ThetaP; Update.VisVec = VisVecP;
end