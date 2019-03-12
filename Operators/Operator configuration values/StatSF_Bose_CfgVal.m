% --- Density correlation configuration value function ---

function [CfgVal] = StatSF_Bose_CfgVal(Hilbert,Cfg,~,~,~)
% Given a bosonic configuration Cfg this function computes the value of the
% static structure factor Nq = <Cfg|n{q}n{-q}|Cfg>. Generalisation to more
% than one dimension may require adjustment.
[Cfg_vec] = Hilbert.FullCfgRef(Cfg);
N = Cfg.N;
CfgVal = zeros(N,1);
CfgQ_vec = fft(Cfg_vec);
for q = 1:N
    CfgVal(q) = CfgQ_vec(1+mod(q-1,N)) * CfgQ_vec(1+mod(1-q,N));
end