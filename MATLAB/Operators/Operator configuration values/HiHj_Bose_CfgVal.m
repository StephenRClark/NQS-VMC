% --- Two site density correlation configuration value function ---

function [CfgVal] = HiHj_Bose_CfgVal(Hilbert,Cfg,~,~,Bonds)
% Given a bosonic configuration Cfg this function computes the matrix
% element <Cfg|[D{i}D{j}]|Cfg>.

[Cfg_vec] = Hilbert.FullCfg(Cfg); % Build the configuration vector for Cfg for convenience below.
N = Cfg.N; % Number of sites in the system.
CoOrd = size(Bonds,2); % Coordination of the lookup list provided.

% This operator is diagonal in the boson number basis.
CfgVal = 0;

for c = 1:CoOrd
    for n=1:N
        m = Bonds(n,c); % Nearest-neighbour site to m.
        CfgVal = CfgVal + ((Cfg_vec(n)==0) * (Cfg_vec(m)==0));
    end
end