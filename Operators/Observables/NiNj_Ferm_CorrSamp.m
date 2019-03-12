% --- Two site density correlation ---

function NiNjSamp = NiNj_Ferm_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the estimator of operator N{i}*N{j} for a
% fermionic configuration. Structure differs from other two site 
% correlators as operator is diagonal in basis.

N = Cfg.N; NiNjSamp = zeros(N,2);

Cfg_vec = (BiPtFermCfg(Cfg) + ones(2*N,1))/2;

for n = 1:N
    for m = 1:N
        NiNjSamp(n,1) = NiNjSamp(n,1) + Cfg_vec(m) * Cfg_vec(1+mod(m+n-2,N));
        NiNjSamp(n,2) = NiNjSamp(n,2) + Cfg_vec(m+N) * Cfg_vec(1+mod(m+n-2,N)+N);
    end
end

NiNjSamp = NiNjSamp / N;