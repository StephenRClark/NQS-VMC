% --- Two-site spin correlation ---

function SiSjSamp = SpSm_Ferm_1D_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the local estimator of operator S+{i}S-{j}  

MCMC.SCorr = @SpSm_Ferm_1D_PBC_CorrMatEls;

SiSjSamp = SpinCorr2S_1D_Samp(Cfg,EnLoc,dLogp,MCMC,Ansatz);