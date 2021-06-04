% --- Two site spin correlation configuration value function ---

function [CfgVal] = SzSz_CfgVal(Hilbert,Cfg,~,~,Bonds)
% Given a spin-1/2 configuration Cfg this function computes the list of
% configurations CfgP's and matrix elements <Cfg|Sz{i}*Sz{j}|CfgP>.

[Cfg_vec] = Hilbert.FullCfg(Cfg); % Build the configuration vector for Cfg for convenience below.
N = Cfg.N; % Number of sites in the system.
CoOrd = size(Bonds,2); % Coordination of the lookup list provided.

% This operator is diagonal in the Sz basis.
CfgVal = 0;

for c = 1:CoOrd
    for n=1:N
        m = Bonds(n,c); % Nearest-neighbour site to m.
        if m ~= 0
            CfgVal = CfgVal + (1/4)*Cfg_vec(n)*Cfg_vec(m); % Compute matrix element.
        end
    end
end