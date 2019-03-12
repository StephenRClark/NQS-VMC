% --- Single site correlation ---

function MSqSamp = LocMoment_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the local moment correlator of the Fermi Hubbard model. 

Cfg_vec = FullFermCfg(Cfg);

MSqSamp = sum(mod(Cfg_vec(:),2))/MCMC.N;