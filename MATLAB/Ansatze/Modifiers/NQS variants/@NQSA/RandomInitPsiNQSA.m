% --- General NQS wave function random initialisation function ---

function [NQSObj] = RandomInitPsiNQSA(NQSObj,Params)
% This function populates random initial NQS ansatz structure. The input
% NQS is assumed to have Nv and Nh defined already. The Params structure
% contains information controlling the form of random elements generated.
% ---------------------------------
% Format for NQSA Modifier:
% - NQSA.Nv = number of "visible" spins.
% - NQSA.Nh = number of "hidden" spins.
% - NQSA.Np = number of parameters in the ansatz = Alpha + Alpha*Nv + 2*Nsl.
% - NQSA.a = (Nv x 1) vector - visible site bias.
% - NQSA.av = (Nsl x 1) vector - visible bias parameters.
% - NQSA.A = (Nv x 1) vector - visible site square bias.
% - NQSA.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSA.b = (Nh x 1) vector - hidden site bias.
% - NQSA.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSA.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQSA.Wm = (Alpha x Nv) matrix - coupling parameters.
% - NQSA.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSA.Theta = (Nh x 1) vector - effective angles.
% - NQSA.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
% Extract information on translational symmetries from Graph.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; SLInds = GraphObj.SLInds;
Nsl = max(SLInds); Ntr = numel(BondMap); Ng = GraphObj.N; Alpha = NQSObj.Alpha;

NQSObj.Np = 2*Nsl + Alpha + (Nv * Alpha); % The number of variational parameters.
Nh = Alpha * Ntr; NQSObj.Nh = Nh;

% Initialise the storage:
NQSObj.a = zeros(Nv,1);
NQSObj.av = zeros(Nsl,1);
NQSObj.A = zeros(Nv,1);
NQSObj.Av = zeros(Nsl,1);
NQSObj.b = zeros(Nh,1);
NQSObj.bv = zeros(Alpha,1);
NQSObj.B = zeros(Nh,1);
NQSObj.Bv = zeros(Alpha,1);
NQSObj.W = zeros(Nh,Nv);
NQSObj.Wm = zeros(Alpha,Nv);
NQSObj.Theta = zeros(Nh,1);

if isfield(Params,'A') == 0
    Params.A = Params.a;
end
if isfield(Params,'B') == 0
    Params.B = Params.b;
end

for v = 1:Nsl
    NQSObj.av(v) = (Params.a + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.a~=0);
    NQSObj.Av(v) = (Params.A + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.A~=0);
end
for al = 1:Alpha
    NQSObj.bv(al) = (Params.b + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.b~=0);
    NQSObj.Bv(al) = (Params.B + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.B~=0);
    for v = 1:Nv
        NQSObj.Wm(al,v) = (Params.W + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.W~=0);
    end
end

% Repackage the ati, bti and Wv to usual NQS form.
for n = 1:Nv
NQSObj.a(n) = NQSObj.av(SLInds(n));
NQSObj.A(n) = NQSObj.Av(SLInds(n));
end
% Constructing shift invariant W matrix.
for al = 1:Alpha
    NQSObj.b((1:Ntr)+(al-1)*Ntr) = NQSObj.bv(al);
    NQSObj.B((1:Ntr)+(al-1)*Ntr) = NQSObj.Bv(al);
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

NQSObj.OptInds = [(real(NQSObj.av)~=0), (imag(NQSObj.av)~=0); (real(NQSObj.Av)~=0),...
    (imag(NQSObj.Av)~=0); (real(NQSObj.bv)~=0), (imag(NQSObj.bv)~=0);...
    (real(NQSObj.Wm(:))~=0), (imag(NQSObj.Wm(:))~=0)];

end