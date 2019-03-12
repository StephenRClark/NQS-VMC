% --- Single site spin correlation matrix element function ---

function [Diff, OpMatEls] = Sm_OpMatEls(HilbertObj,Cfg)
% Given a spin-1/2 configuration Cfg this function computes the list of
% configurations CfgP's and matrix elements <Cfg|S-{i}|CfgP>. Only the
% limited number of "reachable" configurations, where the matrix element is
% non-zero, are returned.

[Cfg_vec] = HilbertObj.FullCfgRef(Cfg); % Build the configuration vector for Cfg for convenience below.
N = Cfg.N; % Number of sites in the system.

% There will be at most N reachable configurations with PBC.
OpMatEls = zeros(N,1);
% Initialise a structure array by populating the final element first:
Diff(N).pos = 1; % To be safe put one site in..
Diff(N).val = 0; %
Diff(N).num = 0; % No differences.
Diff(N).sign = 1; % For compatibility with fermionic implementations.
Diff(N).type = 0; % Also for compatibility with fermionic implementations.
% Now initialise the rest:
for n=1:N
    Diff(n).pos = 0; % Array of positions of configuration differences.
    Diff(n).val = 0; % Array of difference values at those positions.
    Diff(n).num = 1; % Spin ladder operators create only one difference.
    Diff(n).sign = 1;
    Diff(n).type = 0;
end

% Deal with each Sm contribution:
% Compute all configurations obtained from Cfg by a single ladder operator
% Loop over sites and note the difference to the configuration:
for n=1:N
    % Compute the differences CfgP - Cfg:
    
    Diff(n).pos = n;
    Diff(n).val = -2*Cfg_vec(n); % Flip the spin at site (i,j).
    OpMatEls(n) = (1/2) * (Cfg_vec(n)==1); % Compute matrix element
end

% Remove any zero matrix element contributions from the list:
ind = OpMatEls~=0;
OpMatEls = OpMatEls(ind);
Diff = Diff(ind); % Trim the structure array down.