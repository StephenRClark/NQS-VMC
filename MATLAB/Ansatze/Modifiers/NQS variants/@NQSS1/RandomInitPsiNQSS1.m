% --- General NQS wave function random initialisation function ---

function [NQSObj] = RandomInitPsiNQSS1(NQSObj,Params)
% This function populates random initial NQS ansatz structure. The input
% NQS is assumed to have Nv and Nh defined already. The Params structure
% contains information controlling the form of random elements generated.
% ---------------------------------
% Format for NQS Modifier object modified for spin-1:
% - NQSS1.Nv = number of "visible" spins.
% - NQSS1.Nh = number of "hidden" spins.
% - NQSS1.Alpha = number of unique coupling sets or "hidden unit density"
% - NQSS1.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQSS1.a = (Nv x 1) vector - visible site bias.
% - NQSS1.av = (Nsl x 1) vector - visible bias parameters.
% - NQSS1.A = (Nv x 1) vector - visible site square bias.
% - NQSS1.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSS1.b = (Nh x 1) vector - hidden site bias.
% - NQSS1.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSS1.w = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQSS1.wm = (Alpha x Nv) matrix - linear coupling parameters.
% - NQSS1.W = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQSS1.Wm = (Alpha x Nv) matrix - square coupling parameters.
% - NQSS1.Theta = (Nh x 1) vector - effective angles.
% - NQSS1.VisVec = (Nv x 1) vector - visible occupancies vector.
% - NQSS1.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
    
% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" units.
Alpha = NQSObj.Alpha; % Density of "hidden" units.

% Extract information on translational symmetries from Graph.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; Ng = GraphObj.N; SLInds = GraphObj.SLInds;
Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
Nsl = max(SLInds); % Number of sublattices for da.

NQSObj.Np = 2*Nsl + Alpha + 2*(Alpha * Nv); % The number of variational parameters.
Nh = Alpha * Ntr; NQSObj.Nh = Nh; % The number of hidden units, including translates.

% Initialise the storage:
NQSObj.a = zeros(Nv,1);
NQSObj.av = zeros(Nsl,1);
NQSObj.A = zeros(Nv,1);
NQSObj.Av = zeros(Nsl,1);
NQSObj.b = zeros(Nh,1);
NQSObj.bv = zeros(Alpha,1);
NQSObj.w = zeros(Nh,Nv);
NQSObj.wm = zeros(Alpha,Nv);
NQSObj.W = zeros(Nh,Nv);
NQSObj.Wm = zeros(Alpha,Nv);
NQSObj.Theta = zeros(Nh,1);

if isfield(Params,'A') == 0
    Params.A = Params.a;
end
if isfield(Params,'w') == 0
    Params.w = Params.W;
end

for v = 1:Nsl
    NQSObj.av(v) = (Params.a + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.a~=0);
    NQSObj.Av(v) = (Params.A + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.A~=0);
end
for al=1:Alpha
    NQSObj.bv(al) = (Params.b + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.b~=0);
    for v = 1:Nv
        NQSObj.wm(al,v) = (Params.w + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.w~=0);
        NQSObj.Wm(al,v) = (Params.W + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.W~=0);
    end
end

% Repackage the ati, bti and Wv to usual NQS form.
for n = 1:Nv
    NQSObj.a(n) = NQSObj.av(SLInds(n));
    NQSObj.A(n) = NQSObj.Av(SLInds(n));
end
% Constructing shift invariant W matrix.
for a = 1:Alpha
    NQSObj.b((1:Ntr)+(a-1)*Ntr) = NQSObj.bv(a);
    % For each layer labelled by a, generate the desired translates.
    for b = 1:numel(BondMap)
        for n = 1:Nv
            if BondMap{b}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{b}(1+mod(n-1,Ng)) + Ng*(ceil(n/Ng)-1);
                % Account for enlarged lattices where Nv = Ns x Ng.                
                NQSObj.w(b+(a-1)*Ntr,VInd) = NQSObj.wm(a,n);
                NQSObj.W(b+(a-1)*Ntr,VInd) = NQSObj.Wm(a,n);
            end
        end
    end
end

% Set optimisation indicies.
OptInds_w_R = real(NQSObj.wm.')~=0; OptInds_w_I = imag(NQSObj.wm.')~=0;
OptInds_W_R = real(NQSObj.Wm.')~=0; OptInds_W_I = imag(NQSObj.Wm.')~=0;

NQSObj.OptInds = [(real(NQSObj.av)~=0), (imag(NQSObj.av)~=0); (real(NQSObj.Av)~=0),...
    (imag(NQSObj.Av)~=0); (real(NQSObj.bv)~=0), (imag(NQSObj.bv)~=0);...
    OptInds_w_R, OptInds_w_I; OptInds_W_R, OptInds_W_I];

end