% --- Single site spin correlation configuration value function ---

function [CfgVal] = Sz_Sector_CfgVal(Hilbert,Cfg,~,~,~)
% Given a spin-1/2 configuration Cfg this function computes the total Sz
% projection and records it. Output is a vector of Sz sector occupation
% probabilities.

[Cfg_vec] = Hilbert.FullCfg(Cfg); % Build the configuration vector for Cfg for convenience below.

% Sz total can take values: -N, -N+2, ... , N-2, N.

% Total number of possible sectors: N + mod(1+N,2).

SzT = sum(Cfg_vec);

CfgVal = zeros(Cfg.N+mod(Cfg.N+1,2),1); CfgVal(1+(SzT+Cfg.N)/2) = 1;
end