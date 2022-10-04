% --- General NQS wave function random initialisation function ---

function [NQSObj] = RandomInitPsiNQSM(NQSObj,Params)
% This function populates random initial NQS ansatz structure. The input
% NQS is assumed to have Nv and Nh defined already. The Params structure
% contains information controlling the form of random elements generated.
% ---------------------------------
% Format for NQSM Modifier object:
% - NQSM.Nv = number of "visible" units.
% - NQSM.Nh = number of "hidden" units.
% - NQSM.Np = number of parameters in the ansatz = 3*Nv + Alpha + (2*Nv * Alpha).
% - NQSM.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSM.VDim = dimensions of the visible units.
% - NQSM.a = (3*Nv x 1) vector - visible site bias.
% - NQSM.av = (3*Nsl x 1) vector - visible bias parameters.
% - NQSM.b = (Nh x 1) vector - hidden site bias.
% - NQSM.bv =  (Alpha x 1) vector - hidden bias parameters.
% - NQSM.W = (Nh x Nv) matrix - holon coupling terms.
% - NQSM.Wm = (Alpha x Nv) matrix - holon coupling parameters.
% - NQSM.X = (Nh x Nv) matrix - doublon coupling terms.
% - NQSM.Xm = (Alpha x Nv) matrix - doublon coupling parameters.
% - NQSM.Theta = (Nh x 1) vector - effective angles.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; Nmax = NQSObj.VDim-1; % Number and dimension of "visible" units.
Alpha = NQSObj.Alpha; % Density of "hidden" units.

% Extract information on translational symmetries from Graph.
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap; Ng = GraphObj.N; SLInds = GraphObj.SLInds;
Ntr = numel(BondMap); % Number of translates - Nh = Ntr*Alpha.
Nsl = max(SLInds); % Number of sublattices for da.

NQSObj.Np = 3*Nsl + Alpha + (2*Nv * Alpha); % The number of variational parameters.
Nh = Alpha * Ntr; NQSObj.Nh = Nh; % The number of hidden units, including translates.

if ~isfield(Params,"X")
    Params.X = Params.W;
end

% Initialise the storage:
NQSObj.av = zeros(3*Nsl,1);
NQSObj.a = zeros(3*Nv,1);
NQSObj.bv = zeros(Alpha,1);
NQSObj.b = zeros(Nh,1);
NQSObj.Wm = zeros(Alpha,Nv);
NQSObj.W = zeros(Nh,Nv);
NQSObj.Xm = zeros(Alpha,Nv);
NQSObj.X = zeros(Nh,Nv);
NQSObj.Theta = zeros(Nh,1);

for v = 1:3*Nsl
    NQSObj.av(v) = (Params.a + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.a~=0);
end
for h=1:Alpha
    NQSObj.bv(h) = (Params.b + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.b~=0);
end
for h = 1:Alpha
    for v = 1:Nv
        NQSObj.Wm(h,v) = (Params.W + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.W~=0);
        NQSObj.Xm(h,v) = (Params.W + 2*Params.nmag*(rand-0.5)) * exp(2i*pi*Params.nphs*rand)*(Params.X~=0);
    end
end

% Repackage the parameters in the necessary form.
a_h = zeros(Nv,1); a_d = zeros(Nv,1); a_m = zeros(Nv,1);
for s = 1:Nsl
    a_h(SLInds==s) = NQSObj.av(s);
    a_d(SLInds==s) = NQSObj.av(s+Nsl);
    a_m(SLInds==s) = NQSObj.av(s+2*Nsl);
end
NQSObj.a = [a_h; a_d; a_m];
for al = 1:Alpha
    NQSObj.b((1:Ntr)+(al-1)*Ntr) = NQSObj.bv(al);
    % For each layer labelled by a, generate the desired translates.
    for bd = 1:Ntr
        for n = 1:Nv
            if BondMap{bd}(1+mod(n-1,Ng)) ~= 0 % Check that bond is valid - W(b,n) left empty otherwise.
                VInd = BondMap{bd}(1+mod(n-1,Ng));
                % Account for enlarged lattices where Nv = Ns x Ng.
                NQSObj.W(bd+(al-1)*Ntr,VInd) = NQSObj.Wm(al,n);
                NQSObj.X(bd+(al-1)*Ntr,VInd) = NQSObj.Xm(al,n);
            end
        end
    end
end

% Set optimisation indicies.
OptInds_W_R = real(NQSObj.Wm.')~=0; OptInds_W_I = imag(NQSObj.Wm.')~=0;
OptInds_X_R = real(NQSObj.Xm.')~=0; OptInds_X_I = imag(NQSObj.Xm.')~=0;
NQSObj.OptInds = [(real(NQSObj.av)~=0), (imag(NQSObj.av)~=0); (real(NQSObj.bv)~=0),...
    (imag(NQSObj.bv)~=0); OptInds_W_R(:), OptInds_W_I(:); OptInds_X_R(:), OptInds_X_I(:)];

end