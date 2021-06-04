% --- Two site density correlation configuration value function ---

function [CfgVal] = NiNj_Ferm_CfgVal(Hilbert,Cfg,~,~,Bonds)
% Given a fermionic configuration Cfg this function computes the matrix
% element <Cfg|N{i}*N{j}|Cfg>.

[Cfg_vec] = Hilbert.FullCfg(Cfg); % Build the configuration vector for Cfg for convenience below.
N = Cfg.N; % Number of sites in the system.
CoOrd = size(Bonds,2); % Coordination of the lookup list provided.

% This operator is diagonal in the fermionic occupation basis.
CfgVal = 0;

for c = 1:CoOrd
    for n=1:N
        m = Bonds(n,c); % Nearest-neighbour site to m.
        if m ~= 0
            % Correlations covered here are purely density-density.
            CfgVal = CfgVal + sum(Cfg_vec(n,:))*sum(Cfg_vec(m,:)); % Compute matrix element.
        end
    end
end