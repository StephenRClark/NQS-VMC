% --- Two site density correlation ---

function NiNjSamp = NiNj_Bose_CorrSamp(Cfg,EnLoc,dLogp,MCMC,Ansatz)
% This function evaluates the estimator of operator N{i}*N{j} for a
% bosonic configuration. Structure differs from other two site 
% correlators as operator is diagonal in basis.

N = Cfg.N; NiNjSamp = zeros(1,N);

Cfg_vec = Cfg.occ;

for n = 1:N
    for m = 1:N
        NiNjSamp(n) = NiNjSamp(n) + Cfg_vec(m) * Cfg_vec(1+mod(m+n-2,N));
    end
end

NiNjSamp = NiNjSamp / N;