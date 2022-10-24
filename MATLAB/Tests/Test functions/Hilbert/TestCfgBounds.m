function [lower_cond,upper_cond,local_cond,bad_cfgs] = TestCfgBounds(Hilbert,N_cfgs,LowerBound,UpperBound)
% Test function to ensure the configurations generated conform to the local
% and configuration wide restrictions of the Hilbert space. If any are not
% met, then the configurations that fail are outputted.
N = Hilbert.N; d = Hilbert.d; bad_cfgs = zeros(N_cfgs,N); good_cfg_inds = zeros(N_cfgs,1);
lower_cond = true; upper_cond = true; local_cond = true;
n_fail_lbc = 0; n_fail_ubc = 0; n_fail_sbc = 0;
% Generate random configurations:
for n = N_cfgs:1
    Cfg = Hilbert.RandomCfg();
    Cfg_vec = Hilbert.FullCfg(Cfg);
    Cfg_sum = sum(Cfg_vec);
    % Test total particle number / spin are within expected bounds.
    lb_cfg_cond = (Cfg_sum>=LowerBound);
    ub_cfg_cond = (Cfg_sum<=UpperBound);
    if ~lb_cfg_cond
        n_fail_lbc = n_fail_lbc + 1;
    end
    if ~ub_cfg_cond
        n_fail_ubc = n_fail_ubc + 1;
    end
    if strcmp(Hilbert.Type,'Bose')
        lb_site = 0; ub_site = d - 1;
    elseif strcmp(Hilbert.Type,'Spin')
        % NB: Half-integer spin give integer configurations.
        lb_site = -(d-1)/(1+mod(d,2)); ub_site = -lb_site;
    end
    site_cond = true; % Test local values are within expected bounds.
    for s = 1:N
        site_cond = site_cond && (Cfg_vec>=lb_site) && (Cfg_vec<=ub_site);
    end
    if ~site_cond
        n_fail_sbc = n_fail_sbc + 1;
    end
    lower_cond = lower_cond && lb_cfg_cond;
    upper_cond = upper_cond && ub_cfg_cond;
    local_cond = local_cond && site_cond;
    if ~(lb_cfg_cond && ub_cfg_cond && site_cond)
        bad_cfgs(n) = Cfg_vec.';
    else
        good_cfg_inds(n) = 1;
    end
end
if ~lower_cond
    disp(['Proportion of configurations breaching lower total bound: ' ...
        num2str(n_fail_lbc) ' / ' num2str(N_cfgs)]);
end
if ~upper_cond
    disp(['Proportion of configurations breaching upper total bound: ' ...
        num2str(n_fail_ubc) ' / ' num2str(N_cfgs)]);
end
if ~local_cond
    disp(['Proportion of configurations breaching local bounds: ' ...
        num2str(n_fail_sbc) ' / ' num2str(N_cfgs)]);
end
% Trim good cfgs from list.
bad_cfgs = bad_cfgs(good_cfg_inds==0,:);
end