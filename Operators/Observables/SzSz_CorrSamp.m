% --- Two site spin correlation ---

function SiSjSamp = SzSz_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the estimator of operator Sz{i}*Sz{j} - structure
% differs from other two site correlators as operator is diagonal in basis.

N = Cfg.N; SiSjSamp = zeros(N,1);

Cfg_vec = FullSpinCfg(Cfg);

for n = 1:N
    for m = 1:N
        SiSjSamp(n) = SiSjSamp(n) + Cfg_vec(m) * Cfg_vec(1+mod(m+n-2,N));
    end
end

SiSjSamp = SiSjSamp / N;