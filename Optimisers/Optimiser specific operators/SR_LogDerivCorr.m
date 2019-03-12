% --- Stochastic Reconfiguration log-derivative correlation ---

function dLogpCorr = SR_LogDerivCorr(Operator,Cfg,EnLoc,dLogp,Ansatz)
% This function evaluates the log-derivative correlations with itself required
% for Stochastic Reconfiguration.  

dLogpCorr = conj(dLogp) * dLogp.';