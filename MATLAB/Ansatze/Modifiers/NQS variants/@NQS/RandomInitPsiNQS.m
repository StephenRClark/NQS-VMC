% --- General NQS wave function random initialisation function ---

function [NQSObj] = RandomInitPsiNQS(NQSObj,Params)
% This function populates random initial NQS ansatz structure. The input
% NQS is assumed to have Nv and Nh defined already. The Params structure
% contains information controlling the form of random elements generated.
% ---------------------------------
% Format for NQS Modifier object:
% - NQS.Nv = number of "visible" units.
% - NQS.Nh = number of "hidden" units.
% - NQS.Np = number of parameters in the ansatz = Nv + Nh + (Nv * Nh).
% - NQS.Alpha = number of unique coupling sets or "hidden unit density".
% - NQS.a = (Nv x 1) vector - visible site bias.
% - NQS.av = (Nsl x 1) vector - visible bias parameters.
% - NQS.b = (Nh x 1) vector - hidden site bias.
% - NQS.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQS.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQS.Wm = (Alpha x Nv) matrix - hidden-visible coupling parameters.
% - NQS.Theta = (Nh x 1) vector - effective angles.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" units.
Alpha = NQSObj.Alpha; % Density of "hidden" units.

% Extract information on translational symmetries from Graph.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; Ng = GraphObj.N; SLInds = GraphObj.SLInds;
Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
Nsl = max(SLInds); % Number of sublattices for da.

NQSObj.Np = Nsl + Alpha + (Alpha * Nv); % The number of variational parameters.
Nh = Alpha * Ntr; NQSObj.Nh = Nh; % The number of hidden units, including translates.

% Initialise the storage:
NQSObj.a = zeros(Nv,1);
NQSObj.av = zeros(Nsl,1);
NQSObj.b = zeros(Nh,1);
NQSObj.bv = zeros(Alpha,1);
NQSObj.W = zeros(Nh,Nv);
NQSObj.Wm = zeros(Alpha,Nv);
NQSObj.Theta = zeros(Nh,1);

for s = 1:Nsl
    NQSObj.av(s) = (Params.a + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.a~=0);
end
for al=1:Alpha
    NQSObj.bv(al) = (Params.b  + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.b~=0);
end
for al = 1:Alpha
    for v = 1:Nv
        NQSObj.Wm(al,v) = (Params.W  + 2 * Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.W~=0);
    end
end

% Repackage the ati, bti and Wv to usual NQS form.
for n = 1:Nv
    NQSObj.a(n) = NQSObj.av(SLInds(n));
end
% Constructing shift invariant W matrix.
for al = 1:Alpha
    NQSObj.b((1:Ntr)+(al-1)*Ntr) = NQSObj.bv(al);
    % For each layer labelled by a, generate the desired translates.
    for b = 1:numel(BondMap)
        for n = 1:Nv
            if BondMap{b}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{b}(1+mod(n-1,Ng)) + Ng*(ceil(n/Ng)-1);
                % Account for enlarged lattices where Nv = Ns x Ng.
                NQSObj.W(b+(al-1)*Ntr,VInd) = NQSObj.Wm(al,n);
            end
        end
    end
end

% Set optimisation indicies.
OptInds_W_R = real(NQSObj.Wm.')~=0; OptInds_W_I = imag(NQSObj.Wm.')~=0;
NQSObj.OptInds = [(real(NQSObj.av)~=0), (imag(NQSObj.av)~=0); ...
    (real(NQSObj.bv)~=0), (imag(NQSObj.bv)~=0); OptInds_W_R(:), OptInds_W_I(:)];

end