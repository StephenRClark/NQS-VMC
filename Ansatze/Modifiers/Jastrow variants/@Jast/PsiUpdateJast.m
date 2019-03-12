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
% N.B. In fermionic case, density-density terms now have spin indices,
% effectively doubling the size of the lattice (N --> 2N).
% ---------------------------------
% N.B: Jastrow factors here are assumed symmetric i.e.
% Js(i,j) = JS(j,i).
% ---------------------------------
% Format for dLogp vector is a (Np x 1) vector of relevant two-site terms.
% ---------------------------------

JastObj.JsVar = JastObj.JsVar + dP;

cap = JastObj.ParamCap;

% Sanity check the values of the ansatz:
JastObj.JsVar(isinf(JastObj.JsVar)) = 0;
JastObj.JsVar(isnan(JastObj.JsVar)) = 0;
ind = abs(real(JastObj.JsVar))>cap;
JastObj.JsVar(ind) = sign(real(JastObj.JsVar(ind)))*cap + 1i*imag(JastObj.JsVar(ind));

% Repopulate Jast.Js with the new JsVar values.
for n = 1:N
    for m = n:N
       JastObj.Js(m,n) = JastObj.JsVar(N*(n-1) - (n*(n-1)/2) + m);
       if m ~= n
           JastObj.Js(n,m) = JastObj.Js(m,n);
       end
    end
end