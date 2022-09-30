% --- General NQS wave function random initialisation function ---

function [NQSObj] = RandomInitPsiNQSP(NQSObj,Params)
% This function populates random initial NQS ansatz structure. The input
% NQS is assumed to have Nv and Nh defined already. The Params structure
% contains information controlling the form of random elements generated.
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

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; VOrder = NQSObj.VOrder; % Number and dimension of "visible" units.
Alpha = NQSObj.Alpha; HOrder = NQSObj.HOrder; % Density and dimension of "hidden" units.

GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; SLInds = GraphObj.SLInds;

Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
Nsl = max(SLInds); % Number of sublattices for da.

NQSObj.Np = Nsl*VOrder + Alpha*HOrder + Alpha*Nv*HOrder*VOrder; % The number of variational parameters.
Nh = Ntr*Alpha; NQSObj.Nh = Nh; % The number of hidden units, including translates.

% Initialise the storage:
NQSObj.a = zeros(Nv,VOrder);
NQSObj.av = zeros(Nsl,VOrder);
NQSObj.b = zeros(Nh,HOrder);
NQSObj.bv = zeros(Alpha,HOrder);
NQSObj.W = zeros(Nh,Nv,HOrder,VOrder);
NQSObj.Wm = zeros(Alpha,Nv,HOrder,VOrder);
NQSObj.Theta = zeros(Nh,HOrder);

for v = 1:Nsl*VOrder
    NQSObj.av(v) = (Params.a + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.a~=0);
end
for al=1:Alpha
    for ho = 1:HOrder
        NQSObj.bv(al,ho) = (Params.b  + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.b~=0);
        for n = 1:Nv
            for vo = 1:VOrder
                NQSObj.Wv(al,n,ho,vo) = (Params.W  + 2 * Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.W~=0);
            end
        end
    end
end

% Set optimisation indicies.

OptInds_a_R = real(NQSObj.av.')~=0; OptInds_a_I = imag(NQSObj.av.')~=0;
OptInds_b_R = real(NQSObj.bv.')~=0; OptInds_b_I = imag(NQSObj.bv.')~=0;
OptInds_W_R = real(permute(NQSObj.Wm,[4 3 2 1]))~=0;
OptInds_W_I = imag(permute(NQSObj.Wm,[4 3 2 1]))~=0;
NQSObj.OptInds = [OptInds_a_R(:), OptInds_a_I(:); OptInds_b_R(:), OptInds_b_I(:);...
    OptInds_W_R(:), OptInds_W_I(:)];

end