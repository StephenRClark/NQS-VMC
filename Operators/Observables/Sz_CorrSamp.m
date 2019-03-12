% --- Single site spin correlation ---

function SiSamp = Sz_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the local estimator of operator Sz 

Cfg_vec = FullSpinCfg(Cfg);

SiSamp = (1/2)*sum(Cfg_vec(:))/MCMC.N;