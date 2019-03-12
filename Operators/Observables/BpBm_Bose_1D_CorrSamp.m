% --- Two-site boson correlation ---

function BiBjSamp = BpBm_Bose_1D_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the local estimator of operator S+{i}S-{j}  

MCMC.SCorr = @BpBm_1D_PBC_CorrMatEls;

BiBjSamp = BoseCorr2S_1D_Samp(Cfg,EnLoc,dLogp,MCMC,Ansatz);