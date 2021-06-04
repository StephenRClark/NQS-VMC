% --- General Jastrow factor wave function random initialisation function ---

function [JastObj] = RandomInitPsiJast(JastObj,Params)
% This function populates a Jastrow Modifier object with random initial
% parameters. The Params structure contains information controlling the
% form of random elements generated.
% ---------------------------------
% Format for Jastrow Modifier object:
% - Jast.N = number of sites (defined on input).
% - Jast.Np = number of variational Jastrow parameters.
% - Jast.Js = (N x N) matrix - field containing all Jastrow factors.
% - Jast.JsVar = (Np x 1) vector - Jastrow variational parameters.
% - Jast.Tj = (N x 1) vector - used to track on-site contributions.
% - Jast.JsV = (N x N) matrix - contains variational parameter indices for each site.
% ---------------------------------
% N.B: Jastrow factors here are assumed symmetric i.e.
% Js(i,j) = JS(j,i).
% ---------------------------------

% Make local copies to reduce notation in code below.
Nj = JastObj.N; % Number of sites for Jastrow factors.
GraphObj = JastObj.Graph; Ng = GraphObj.N; % Number of actual sites in Graph.
BondMap = GraphObj.BondMap;

JsV = zeros(Nj); Np = 0; % Tracker for number of parameters required.
% Translation invariance - Js(i,j) = Js(i+dR,j+dR).
for n = 1:Nj
    for m = n:Nj % Symmetry requirement means initially only need to populate on one side of the diagonal.
        if JsV(n,m) == 0 % If not populated with a variable index, needs assignment.
            Np = Np + 1;
            JsV(n,m) = Np; JsV(m,n) = Np;
            for b = 1:numel(BondMap) % Assign variable indices for each translate.
                Inds = [(BondMap{b}(1+mod(m-1,Ng))+Ng*(ceil(m/Ng)-1)),...
                    (BondMap{b}(1+mod(n-1,Ng))+Ng*(ceil(n/Ng)-1))];
                if Inds(1) ~= 0 && Inds(2) ~= 0
                    JsV(Inds(1),Inds(2)) = Np; JsV(Inds(2),Inds(1)) = Np;
                end
            end
        end
    end
end

% Initialise the storage:
JastObj.Np = Np;
JastObj.Js = zeros(Nj);
JastObj.JsVar = zeros(Np,1);
JastObj.Tj = zeros(Nj,1);
JastObj.JsV = JsV;

for p = 1:JastObj.Np
    JastObj.JsVar(p) = Params.Js * (1 - Params.nmag + 2*Params.nmag*rand);
    % * exp(2i*pi*Params.nphs*rand); % Jastrow factors normally real.
end

JastObj.Js = JastObj.JsVar(JsV); JastObj.OptInds = ones(JastObj.Np,1);