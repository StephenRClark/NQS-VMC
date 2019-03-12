% --- Single site spin correlation ---

function SiSamp = Sx_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the local estimator of operator Sx 

MCMC.SCorr = @Sx_CorrMatEls;

SiSamp = SpinCorr1S_Samp(Cfg,EnLoc,dLogp,MCMC,Ansatz);