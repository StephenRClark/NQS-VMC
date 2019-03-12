% --- Two-site fermion correlation ---

function CiCjSamp = CpCm_Ferm_1D_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the local estimator of operator S+{i}S-{j}  

MCMC.SCorr = {@CpCmUp_1D_PBC_CorrMatEls,@CpCmDn_1D_PBC_CorrMatEls};

CiCjSamp = FermCorr2S_1D_Samp(Cfg,EnLoc,dLogp,MCMC,Ansatz);