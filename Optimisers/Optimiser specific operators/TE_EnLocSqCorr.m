% --- Local energy squared sampler ---

function EnLocSq = TE_EnLocSqCorr(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the square of the local energy - used for
% variance calculations.

EnLocSq = abs(EnLoc)^2;