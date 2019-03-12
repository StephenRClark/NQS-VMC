% --- General Jastrow logarithmic derivative function ---

function dLogp = LogDerivJastTI(JastObj,HilbertObj,Cfg)
% This function computes the logarithmic derivative:
%            dLogp = 1/Psi(Cfg) dPsi(Cfg)/dp
% w.r.t. each parameter p of the Jastrow ansatz, for a bosonic
% configuration specifed by the structure Cfg.
% ---------------------------------
% Format for Jastrow Modifier object:
% - Jast.N = number of sites (defined on input).
% - Jast.Np = number of variational Jastrow parameters.
% - Jast.Js = (N x N) matrix - field containing all Jastrow factors.
% - Jast.JsVar = (Np x 1) vector - Jastrow variational parameters.
% - Jast.Tj = (N x 1) vector - used to track on-site contributions.
% Properties added with translation invariance:
% - Jast.JsV = (N x N) matrix - contains variational parameter indices for each site.
% N.B. In fermionic case, density-density terms now have spin indices,
% effectively doubling the size of the lattice (N --> 2N).
% ---------------------------------
% N.B: Jastrow factors here are assumed symmetric i.e.
% Js(i,j) = JS(j,i).
% ---------------------------------
% Format for dLogp vector is a (Np x 1) vector of relevant two-site terms.
% ---------------------------------

% Construct vector representation of Cfg.
Cfg_vec = sum(HilbertObj.FullCfgRef(Cfg),2);
% Sum over 2nd index accounts for fermionic representations.

dLogp = zeros(JastObj.Np,1); % Initialise full vector of derivatives.

% Parameter indices associated with site pairs contained in JsV.
for p = 1:JastObj.Np
    % Collect all two-site contributions related to parameter p.
    DDMatP = - ((JastObj.JsV == p) .* (Cfg_vec * Cfg_vec'))/2;
    dLogp(p) = sum(DDMatP(:));
end

% Do some forward error prevention for NaN or Inf elements by zeroing them:
dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;