% --- Two-site spin correlation ---

function SiSjSamp = SpSm_FFT_1D_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the local estimator of the finite Fourier 
% transform of operator S+{i}S-{j}  

MCMC.SCorr = @SpSm_1D_PBC_CorrMatEls;

SiSjSamp = SpinCorr2S_FFT_1D_Samp(Cfg,EnLoc,dLogp,MCMC,Ansatz);