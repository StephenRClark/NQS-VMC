% --- Boson addition matrix element function ---

function [Diff, OpMatEls] = BpQ_OpMatEls(HilbertObj,Cfg,~)
% Given a boson configuration Cfg this function computes the list of
% configurations CfgP's and matrix elements
% <Cfg|exp{1i*qj}B+{i}|CfgP>. Only the limited number of "reachable"
% configurations, where the matrix element is non-zero, are returned.

% As operator acts to left, the conjugate is applied to the supplied
% configuration.

[Cfg_vec] = HilbertObj.FullCfg(Cfg); % Build the configuration vector for Cfg for convenience below.
N = Cfg.N; % Number of sites in the system.

% There will be at most N reachable configurations.
OpMatEls = zeros(N); q = ((1:N)-1)*2*pi/N;
% Each new configuration has N associated values (from N values of
% momentum q).

% Initialise a structure array by populating the final element first:
Diff(N).pos = 1; % To be safe put one site in..
Diff(N).val = 0; %
Diff(N).num = 0; % No differences.
Diff(N).type = 0; % For compatibility.
Diff(N).sign = 1; % For compatibility.

% Now initialise the rest:
for n = 1:N
    Diff(n).pos = n;
    Diff(n).val = -1; % Remove boson from site n.
    OpMatEls(n,:) = sqrt(Cfg_vec(n)) * exp(1i*n*q)/sqrt(N);
    % Second Boolean term comes into play for hardcore boson cases.
end

% Remove any zero matrix element contributions from the list:
ind = OpMatEls(:,1)~=0;
OpMatEls = OpMatEls(ind,:);
Diff = Diff(ind); % Trim the structure array down.