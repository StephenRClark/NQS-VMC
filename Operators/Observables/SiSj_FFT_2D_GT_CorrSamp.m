% --- Two-site spin correlation ---

function SiSjSamp = SiSj_FFT_2D_GT_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the local estimator of operator S{i} x S{j}  

MCMC.SCorr = @SiSj_2D_PBC_GT_CorrMatEls;

SiSjSamp = SpinCorr2S_FFT_2D_Samp(Cfg,EnLoc,dLogp,MCMC,Ansatz);