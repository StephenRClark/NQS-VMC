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
% N.B. In fermionic case, density-density terms now have spin indices,
% effectively doubling the size of the lattice (N --> 2N).
% ---------------------------------
% N.B: Jastrow factors here are assumed symmetric i.e.
% Js(i,j) = JS(j,i).
% ---------------------------------

% Make local copies to reduce notation in code below.
N = JastObj.N; % Number of sites.

Np = N*(N+1)/2; % Number of distinct Jastrow factors with symmetry imposed.

% Initialise the storage:
JastObj.Np = Np;
JastObj.Js = zeros(N);
JastObj.JsVar = zeros(Np,1);
JastObj.Tj = zeros(N,1);

for p = 1:JastObj.Np
    JastObj.JsVar(p) = Params.Js * (1 - Params.nmag + 2*Params.nmag*rand);
    % * exp(2i*pi*Params.nphs*rand); % Jastrow factors normally real.
end

for n = 1:N
    for m = n:N
       JastObj.Js(m,n) = JastObj.JsVar(N*(n-1) - (n*(n-1)/2) + m);
       if m ~= n
           JastObj.Js(n,m) = JastObj.Js(m,n);
       end
    end
end