% --- Occupation number fraction configuration value function ---

function [CfgVal] = OccFrac_Bose_CfgVal(Hilbert,Cfg,~,~,~)
% Given a bosonic configuration Cfg, this function outputs a vector of
% occupation number fractions (fraction of empty, 1-boson, 2-boson etc.
% sites) up to a truncation provided by HilbertObj.

CfgVal = zeros(Hilbert.d,1); N = Hilbert.N;

Cfg_vec = Hilbert.FullCfg(Cfg);

for n = 0:(Hilbert.d-1)
    CfgVal(n+1) = sum(Cfg_vec==n)/N;
end

end