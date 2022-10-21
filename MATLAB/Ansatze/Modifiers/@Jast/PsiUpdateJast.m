% --- General Jastrow wave function update function ---

function JastObj = PsiUpdateJast(JastObj,dP)
% This function updates the Jastrow parameters of the ansatz from a vector
% of parameters dP.
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
% Format for dLogp vector is a (Np x 1) vector of relevant two-site terms.
% ---------------------------------

dP = real(dP).*JastObj.OptInds(:,1) + 1i*imag(dP).*JastObj.OptInds(:,2); % Zeroes out any undesired parameter changes.

JastObj.JsVar = JastObj.JsVar + dP;

cap = JastObj.ParamCap;

% Sanity check the values of the ansatz:
JastObj.JsVar(isinf(JastObj.JsVar)) = 0;
JastObj.JsVar(isnan(JastObj.JsVar)) = 0;
ind = abs(real(JastObj.JsVar))>cap;
JastObj.JsVar(ind) = sign(real(JastObj.JsVar(ind)))*cap + 1i*imag(JastObj.JsVar(ind));

% Normalise so that the sum of Js is zero if flag is active.
if JastObj.NormFlag == 1
    JastObj.JsVar = JastObj.JsVar - mean(JastObj.JsVar(JastObj.JsV(:)));
end
% Repopulate Jast.Js with the new JsVar values.
JastObj.Js = JastObj.JsVar(JastObj.JsV);
end