% --- Two-site spin correlation ---

function SiSjSamp = SpSm_2D_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the local estimator of operator S+{i}S-{j}  

MCMC.SCorr = @SpSm_2D_PBC_CorrMatEls;

SiSjSamp = SpinCorr2S_2D_Samp(Cfg,EnLoc,dLogp,MCMC,Ansatz);