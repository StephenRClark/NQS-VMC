% --- Two-site pair hopping correlation ---

function PiPjSamp = PiPj_Ferm_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the estimator of pair hopping operator 
% C*{i,dn}C*{i,up}C{j,up}C{j,dn} for a fermionic configuration. 

MCMC.SCorr = @PpPm_1D_PBC_CorrMatEls;

PiPjSamp = SpinCorr2S_1D_Samp(Cfg,EnLoc,dLogp,MCMC,Ansatz);