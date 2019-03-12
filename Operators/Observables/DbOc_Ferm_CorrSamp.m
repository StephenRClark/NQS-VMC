% --- Single site density-density correlation ---

function DbSamp = DbOc_Ferm_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the estimator of operator N{i,up}N{i,dn} for a
% fermionic configuration. Structure differs from other two site 
% correlators as operator is diagonal in basis.

Cfg_vec = FullFermCfg(Cfg);

DbSamp = sum(Cfg_vec==2)/Cfg.N;