% --- Multi site spin correlation matrix element function ---

function [Diff, OpMatEls] = SMS_XZZX_OpMatEls(HilbertObj,Cfg,GraphObj,ListIndex)
% Given a spin-1/2 configuration Cfg this function computes the list of
% configurations CfgP's and matrix elements <Cfg|Sx{i}Sy{j}Sx{k}|CfgP>. Only
% the limited number of "reachable" configurations, where the matrix
% element is non-zero, are returned.

[Cfg_vec] = HilbertObj.FullCfg(Cfg); % Build the configuration vector for Cfg for convenience below.

List = GraphObj.ExtraLabels{ListIndex}; Ng = size(List,1);

% There will be at most Ngroup reachable configurations.
OpMatEls = zeros(Ng,1);
% Initialise a structure array by populating the final element first:
Diff(Ng).pos = 1; % To be safe put one site in..
Diff(Ng).val = 0; %
Diff(Ng).num = 0; % No differences.
Diff(Ng).sign = 1; % For compatibility with fermionic implementations.
Diff(Ng).type = 0; % Also for compatibility with fermionic implementations.
Diff(Ng).index = Ng; % Group index, for identification in GraphSample.
% Loop through given lists and grab site values from Cfg_vec.
for n=1:Ng
    % Array of positions of configuration differences.
    Diff(n).pos = [List(n,1), List(n,4)];
    % Array of difference values at those positions.
    Diff(n).val = [-2*Cfg_vec(List(n,1)), -2*Cfg_vec(List(n,4))]; % Flip all relevant spins.
    Diff(n).num = 2; % This operator flips each site.
    Diff(n).sign = 1;
    Diff(n).type = 0;
    Diff(n).index = n;
    OpMatEls(n) = Cfg_vec(List(n,2))*Cfg_vec(List(n,3)); % This Operator acts on all configurations.
end

% Remove any zero matrix element contributions from the list:
ind = OpMatEls~=0;
OpMatEls = OpMatEls(ind);
Diff = Diff(ind); % Trim the structure array down.
end