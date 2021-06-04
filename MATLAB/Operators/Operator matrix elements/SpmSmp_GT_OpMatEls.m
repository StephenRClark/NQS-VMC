% --- Two site spin correlation matrix element function ---

function [Diff, OpMatEls] = SpmSmp_GT_OpMatEls(HilbertObj,Cfg,GraphObj)
% Given a configuration Cfg and a set of linked sites in Bonds, this
% function computes the list of configurations CfgP's and matrix elements
% <Cfg|1/2*(S+{i}*S-{j} + S-{i}*S+{j})|CfgP>. Only the limited number of
% "reachable" configurations, where the matrix element is non-zero, are
% returned.

[Cfg_vec] = HilbertObj.FullCfg(Cfg); % Build the configuration vector for Cfg for convenience below.
N = Cfg.N; % Number of sites in the system.
Bonds = GraphObj.Bonds; SLInds = GraphObj.SLInds;
Nbond = numel(Bonds); % Count the number of bonds provided.
CoOrd = size(Bonds,2); % Coordination of the lookup list provided.

OpMatEls = zeros(Nbond,1);
% Initialise a structure array by populating the final element first:
Diff(Nbond).pos = 1; % To be safe put one site in..
Diff(Nbond).val = 0; %
Diff(Nbond).num = 0; % No differences.
Diff(Nbond).sign = 1; % For compatibility with fermionic implementations.
Diff(Nbond).type = 0; % Also for compatibility with fermionic implementations.
% Now initialise the rest:
for n=1:Nbond
    Diff(n).pos = zeros(1,2); % Array of positions of configuration differences.
    Diff(n).val = zeros(1,2); % Array of difference values at those positions.
    Diff(n).num = 2; % Two-spin terms create at most 2 two differences.
    Diff(n).sign = 1;
    Diff(n).type = 0;
end
for c = 1:CoOrd
    for n = 1:N
        m = Bonds(n,c);
        if m == n
            % Operator is equivalent to identity with zero separation...
            Diff(1+(c-1)*N).num = 0;
            Diff(1+(c-1)*N).pos = 1;
            Diff(1+(c-1)*N).val = 0;
            OpMatEls(1+(c-1)*N) = 1;
            % Only one difference required for this set - break to truncate list.
            break
        else
            if Bonds(n) ~= 0
                Diff(n+(c-1)*N).pos(1) = n;
                Diff(n+(c-1)*N).val(1) = -2*Cfg_vec(n); % Flip the spin at site n.
                Diff(n+(c-1)*N).pos(2) = m;
                Diff(n+(c-1)*N).val(2) = -2*Cfg_vec(m); % Flip the spin at site m.
                OpMatEls(n+(c-1)*N) = ((-1)^(SLInds(m) - SLInds(n))) * ...
                    (1/4)*(1 - Cfg_vec(n)*Cfg_vec(m)); % Non-zero only if spins are different.
            end
        end
    end
end

% Remove any zero matrix element contributions from the list:
ind = OpMatEls~=0;
OpMatEls = OpMatEls(ind);
Diff = Diff(ind); % Trim the structure array down.

