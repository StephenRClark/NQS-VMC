% --- Two-site fermion correlation ---

function CiCjSamp = CpCm_Ferm_2D_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the local estimator of operator S+{i}S-{j}  

MCMC.SCorr = {@CpCmUp_2D_PBC_CorrMatEls,@CpCmDn_2D_PBC_CorrMatEls};

CiCjSamp = FermCorr2S_2D_Samp(Cfg,EnLoc,dLogp,MCMC,Ansatz);