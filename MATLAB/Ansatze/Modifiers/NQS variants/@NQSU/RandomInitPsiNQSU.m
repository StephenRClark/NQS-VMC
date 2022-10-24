% --- General NQS wave function random initialisation function ---

function [NQSObj] = RandomInitPsiNQSU(NQSObj,Params)
% This function populates random initial NQS ansatz structure. The input
% NQS is assumed to have Nv and Nh defined already. The Params structure
% contains information controlling the form of random elements generated.
% ---------------------------------
% Format for NQSU Modifier object:
% - NQSU.Nv = number of "visible" units.
% - NQSU.Nh = number of "hidden" units.
% - NQSU.Np = number of parameters in the ansatz = Nmax*Nv + Nh + (Nmax*Nv * Nh).
% - NQSU.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSU.VDim = dimensions of the visible units.
% - NQSU.a = (Nmax*Nv x 1) vector - visible site bias.
% - NQSU.av = (Nmax*Nsl x 1) vector - visible bias parameters.
% - NQSU.b = (Nh x 1) vector - hidden site bias.
% - NQSU.bv =  (Alpha x 1) vector - hidden bias parameters.
% - NQSU.W = (Nh x Nmax*Nv) matrix - hidden-visible coupling terms.
% - NQSU.Wm = (Alpha x Nmax*Nv) matrix - coupling parameters.
% - NQSU.Theta = (Nh x 1) vector - effective angles.
% - NQSU.VList = (VDim x 1) vector - visible site value list for unary encoding.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; Nmax = NQSObj.VDim-1; % Number and dimension of "visible" units.
Alpha = NQSObj.Alpha; % Density of "hidden" units.

% Extract information on translational symmetries from Graph.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; Ng = GraphObj.N; SLInds = GraphObj.SLInds;
Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
Nsl = max(SLInds); % Number of sublattices for da.

NQSObj.Np = Nmax*Nsl + Alpha + (Nmax*Nv * Alpha); % The number of variational parameters.
Nh = Alpha * Ntr; NQSObj.Nh = Nh; % The number of hidden units, including translates.

% Initialise the storage:
NQSObj.av = zeros(Nmax*Nsl,1);
NQSObj.a = zeros(Nmax*Nv,1);
NQSObj.bv = zeros(Alpha,1);
NQSObj.b = zeros(Nh,1);
NQSObj.Wm = zeros(Alpha,Nmax*Nv);
NQSObj.W = zeros(Nh,Nmax*Nv);
NQSObj.Theta = zeros(Nh,1);

for v = 1:Nmax*Nsl
    NQSObj.av(v) = (Params.a + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.a~=0);
end
for h=1:Alpha
    NQSObj.bv(h) = (Params.b + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.b~=0);
end
for h = 1:Alpha
    for v = 1:Nmax*Nv
        NQSObj.Wm(h,v) = (Params.W + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.W~=0);
    end
end

% Repackage the parameters in the necessary form.
for n = 1:Nv
    for v = 1:Nmax
        NQSObj.a(v + (n-1)*Nmax) = NQSObj.av(v+(SLInds(v)-1)*Nmax);
    end
end
for al = 1:Alpha
    NQSObj.b((1:Ntr)+(al-1)*Ntr) = NQSObj.bv(al);
    % For each layer labelled by a, generate the desired translates.
    for bd = 1:Ntr
        for n = 1:Nv
            if BondMap{bd}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{bd}(1+mod(n-1,Ng));
                for v = 1:Nmax
                    NInd = v + Nmax*(n-1); WInd = v + Nmax*(VInd-1);
                    % Account for enlarged lattices where Nv = Ns x Ng.
                    NQSObj.W(bd+(al-1)*Ntr,WInd) = NQSObj.Wm(al,NInd);
                end
            end
        end
    end
end

% Set optimisation indicies.
NQSObj.OptInds = [(real(NQSObj.av)~=0), (imag(NQSObj.av)~=0); (real(NQSObj.bv)~=0),...
    (imag(NQSObj.bv)~=0); (real(NQSObj.Wm(:))~=0), (imag(NQSObj.Wm(:))~=0)];

end