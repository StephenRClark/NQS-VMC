% --- Density fluctuation correlation ---

function VarNSamp = VarN_Bose_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the estimator of density fluctuations 
% <N^2> - <N>^2 for a bosonic configuration. 

N = Cfg.N; VarNSamp = 0;

Cfg_vec = Cfg.occ;

BoseDen = sum(Cfg_vec)/N;

for n = 1:N
    VarNSamp = VarNSamp + (Cfg_vec(n)^2);
end

VarNSamp = (VarNSamp / N) - (BoseDen^2);