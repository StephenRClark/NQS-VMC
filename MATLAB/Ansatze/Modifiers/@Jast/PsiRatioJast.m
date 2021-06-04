% --- General Jastrow wave function amplitude ratio function ---

function [Ratio,TjP] = PsiRatioJast(JastObj,Diff)
% This function computes the ratio Psi(CfgP)/Psi(Cfg) of amplitudes for a
% proposed configuration CfgP and the current on Cfg, whose difference is
% stored in Diff.
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
% Format for Update is a vector of new local Jastrow contributions.
% ---------------------------------

% Make local copies to reduce notation in the code below:
N = JastObj.N; Tj = JastObj.Tj; TjP = Tj;

if Diff.num ~= 0
    % For multiple changes in on-site quantum numbers, local dependency is
    % handled by Jast.Tj terms, and 'leftover' terms are contained in the
    % following:
    
    % Reshape to ensure no issues with different Diff formats.
    DiffVal = reshape(Diff.val,Diff.num,1);
    DiffPos = reshape(Diff.pos,Diff.num,1);
    % Requires the number of elements in DiffVal and DiffPos to be the same -
    % for fermionic treatments, some pre-processing of Diff may be required.
    
    JsRem = JastObj.Js(DiffPos,DiffPos) .* (DiffVal * DiffVal');
    
    Ratio = exp(-sum(JsRem(:))/2);
    for d = 1:Diff.num 
        % Diff.num should reflect number of sites that change local quantum numbers.
        Ratio = Ratio * exp(-Tj(Diff.pos(d)) * Diff.val(d));
        for n = 1:N
            TjP(n) = TjP(n) + Diff.val(d) * JastObj.Js(Diff.pos(d),n);
        end
    end
else
    Ratio = 1;
end