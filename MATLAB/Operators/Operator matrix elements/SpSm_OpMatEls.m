% --- Two site spin correlation matrix element function ---

function [Diff, OpMatEls] = SpSm_OpMatEls(HilbertObj,Cfg,GraphObj)
% Given a configuration Cfg and a set of linked sites in Bonds, this
% function computes the list of configurations CfgP's and matrix elements
% <Cfg|S+{i}*S-{j}|CfgP>. Only the limited number of "reachable"
% configurations, where the matrix element is non-zero, are returned.

% As operator acts to left, the conjugate is applied to the supplied
% configuration.

[Cfg_vec] = HilbertObj.FullCfg(Cfg); % Build the configuration vector for Cfg for convenience below.
N = Cfg.N; % Number of sites in the system.
Bonds = GraphObj.Bonds;
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
            % No differences produced, and value is 1 if spin is in down
            % state - only need one Diff and OpMatEls element.
            Diff(1+(c-1)*N).num = 0;
            Diff(1+(c-1)*N).pos = 1;
            Diff(1+(c-1)*N).val = 0;
            OpMatEls(1+(c-1)*N) = (1/2)*sum(Cfg_vec==-1);
            break
        else
            if m ~= 0
                Diff(n+(c-1)*N).num = 2;
                Diff(n+(c-1)*N).pos(1) = n;
                Diff(n+(c-1)*N).val(1) = -2; % Lower spin at first site.
                Diff(n+(c-1)*N).pos(2) = m;
                Diff(n+(c-1)*N).val(2) = 2; % Raise spin at bonded site.
                OpMatEls(n+(c-1)*N) = (1/2) * (Cfg_vec(n) > 0) * (Cfg_vec(m) < 0);
            end
        end
    end
end
% Remove any zero matrix element contributions from the list:
ind = OpMatEls~=0;
OpMatEls = OpMatEls(ind);
Diff = Diff(ind); % Trim the structure array down.

