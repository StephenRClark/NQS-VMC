% --- Multi site spin correlation matrix element function ---

function [Diff, OpMatEls] = Sz_MS_OpMatEls(HilbertObj,Cfg,GraphObj,ListIndex)
% Given a spin-1/2 configuration Cfg this function computes the list of
% configurations CfgP's and matrix elements <Cfg|Sz{i,...,j}|CfgP>. Only
% the limited number of "reachable" configurations, where the matrix
% element is non-zero, are returned.

% While diagonal, implemented for OperatorMS operators to distinguish
% groups and between ExtraLabel sets.

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
    Diff(n).pos = 1; 
    % Array of difference values at those positions.
    Diff(n).val = 0; 
    Diff(n).num = 0; % This operator flips each site.
    Diff(n).sign = 1;
    Diff(n).type = 0;
    Diff(n).index = n;
    % This Operator acts on all configurations.
    OpMatEls(n) = prod(Cfg_vec(List(n,:))); 
end

% Remove any zero matrix element contributions from the list:
ind = OpMatEls~=0;
OpMatEls = OpMatEls(ind);
Diff = Diff(ind); % Trim the structure array down.
end