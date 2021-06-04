% --- Two site density correlation configuration value function ---

function [CfgVal] = NiNj_Bose_CfgVal(Hilbert,Cfg,~,~,Bonds)
% Given a bosonic configuration Cfg this function computes the matrix
% element <Cfg|N{i}*N{j}|Cfg>.

[Cfg_vec] = Hilbert.FullCfg(Cfg); % Build the configuration vector for Cfg for convenience below.
N = Cfg.N; % Number of sites in the system.
CoOrd = size(Bonds,2); % Coordination of the lookup list provided.

% This operator is diagonal in the boson number basis.
CfgVal = 0;

for c = 1:CoOrd
    for n=1:N
        m = Bonds(n,c); % Nearest-neighbour site to m.
        if m == n % Expectation value changes if on same site.
            CfgVal = CfgVal + Cfg_vec(n)*(Cfg_vec(n)-1); % Compute matrix element.
        elseif m ~= 0
            % Correlations covered here are purely density-density.
            CfgVal = CfgVal + Cfg_vec(n)*Cfg_vec(m); % Compute matrix element.
        end
    end
end