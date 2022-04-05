% --- General NQS wave function amplitude ratio function ---

function [Ratio,Update] = PsiRatioNQSP(NQSObj,Diff)
% This function computes the ratio Psi(CfgP)/Psi(Cfg) of amplitudes for
% a proposed spin-1/2 configuration CfgP and the current on Cfg, whose
% difference is stored in Diff.
% ---------------------------------
% Format for NQS Modifier object with number hidden units:
% - NQS.Nv = number of "visible" spins.
% - NQS.Nh = number of "hidden" spins.
% - NQS.VDim = (1 x 1) scalar - dimension of visible neurons.
% - NQS.HDim = (1 x 1) scalar - dimension of hidden units.
% - NQS.VOrder = (1 x 1) scalar - highest power of visible unit
% interactions. Max value VDim-1.
% - NQS.HOrder = (1 x 1 ) scalar - highest power of hidden unit
% interactions. Max value HDim-1.
% - NQS.Np = number of parameters in the ansatz = (Nv x VOrder) + (Nh x
% HOrder) + (Nv x VOrder)(Nh x HOrder)
% - NQS.a = (Nv x VOrder) matrix - visible site biases.
% - NQS.b = (Nh x HOrder) matrix - hidden site bias.
% - NQS.W = (Nh x Nv x HOrder x VOrder) tensor - hidden-visible coupling terms.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% - NQS.VisVec = (Nv x 1) vector - visible occupancies.
% ---------------------------------
% Format for Update is a struct with two fields:
% Update.Theta - matrix of new effective angles ThetaP.
% Update.VisVec - vector of new rescaled visible occupancies.
% ---------------------------------

VisVecP = NQSObj.VisVec(:); 
VisPow_shift = zeros(Diff.num,NQSObj.VOrder); % Shift in visible occupancy powers for each altered site.
VisPow_shift(:,1) = Diff.val(:)/(NQSObj.VDim-1);
% Calculate power shifts and visible bias contributions.
Ratio = exp(sum(a(Diff.pos(:),1).*VisPow_shift(:,1)));
for vo = 2:NQSObj.VOrder
    VisPow_shift(:,vo) = (VisVecP(Diff.pos(:)+VisPow_shift(:,1))).^vo - VisVecP(Diff.pos(:)).^vo;
    Ratio = Ratio * exp(sum(a(Diff.pos(:),vo).*VisPow_shift(:,vo)));
end
Theta_shift = zeros(NQSObj.Nh,NQSObj.HOrder); % Initialise effective angle shift.
% Theta(h,ho) = sum(W(h,:,ho,vo)*VisVec.^vo)+b(h,ho);
% Only loop over the sites where there are differences:
for d=1:Diff.num
    for ho = 1:NQSObj.HOrder
        Theta_shift(:,ho) = Theta_shift(:,ho) + sum(reshape(...
            NQSObj.W(:,Diff.pos(d),ho,:),NQSObj.Nh,NQSObj.VOrder).*VisPow_shift(d,:),2);
    end
end
ThetaP = NQSObj.Theta + Theta_shift; % Update the effective angles for the proposed configuration.
% Calculate hidden unit contributions to ratio.
HidVals = zeros(NQSObj.HDim,NQSObj.HOrder);
HidVals(:,1) = (0:(NQSObj.HDim-1))/(NQSObj.HDim-1).';
for ho = 2:NQSObj.HOrder
    HidVals(:,ho) = HidVals(:,1).^ho;
end
HidAmp = zeros(NQSObj.Nh,1); HidAmpP = zeros(NQSObj.Nh,1);
for h = 1:NQSObj.HDim
    HidAmp = HidAmp + exp(NQSObj.Theta*(HidVals(h,:).'));
    HidAmpP = HidAmpP + exp(ThetaP*(HidVals(h,:).'));
end
Ratio = Ratio * prod(HidAmpP./HidAmp);
VisVecP(Diff.pos(:)) = VisVecP(Diff.pos(:))+Diff.val(:)/(NQSObj.VDim-1);
% Add new VisVec and Theta to Update.
Update.Theta = ThetaP; Update.VisVec = VisVecP;
end