% --- Stochastic Reconfiguration energy-log-derivative correlation ---

function En_dLogp = SR_EnLogDerivCorr(Operator,Cfg,EnLoc,dLogp,Ansatz)
% This function evaluates the energy-log-derivative correlation required
% for Stochastic Reconfiguration.  

En_dLogp = EnLoc * conj(dLogp);