% --- Single site spin correlation ---

function SzProfSamp = SzProf_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the local Sz spin profile.

SzProfSamp = (1/2)*FullSpinCfg(Cfg);