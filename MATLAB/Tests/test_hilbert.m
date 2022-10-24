addpath(genpath("../Hilbert")); addpath("Test functions/Hilbert/");

% Set some basic test cases:

Spin_half_free = Spin(10,1/2,[]); % Spin 1/2, no total Sz restriction.
Spin_half_zero = Spin(10,1/2,0); % Spin 1/2, total Sz fixed to zero.
Spin_one_free = Spin(10,1,[]); % Spin 1, no total Sz restriction.
Spin_one_zero = Spin(10,1,0); % Spin 1, total Sz fixed to zero.
Bose_hardcore_free = Bose(10,[],1); % Bosons, max occupation 1, no number restriction.
Bose_hardcore_half = Bose(10,5,1); % Bosons, max occupation 1, half filling.
Bose_4max_free = Bose(10,[],4); % Bosons, max occupation 4, no number restriction.
Bose_4max_den1 = Bose(10,10,4); % Bosons, max occupation 4, average density 1.
Bose_4max_den2 = Bose(10,20,4); % Bosons, max occupation 4, average density 2.

HilbertArray = {Spin_half_free; Spin_half_zero; Spin_one_free; Spin_one_zero; ...
    Bose_hardcore_free; Bose_hardcore_half; Bose_4max_free; Bose_4max_den1; Bose_4max_den2};

HilbertIDs = {'Spin 1/2, free total Sz'; 'Spin 1/2, zero total Sz'; 'Spin 1, free total Sz'; ...
    'Spin 1, zero total Sz'; 'Bosons, hardcore, free total Nb'; 'Bosons, hardcore, half filled'; ... 
    'Bosons, max 4 per site, free total Nb'; 'Bosons, max 4 per site, fixed density 1';...
    'Bosons, max 4 per site, fixed density 2'};

Cfg_lower_bounds = [-10; 0; -10; 0; 0; 5; 0; 10; 20];
Cfg_upper_bounds = [10; 0; 10; 0; 10; 5; 40; 10; 20];

N_cfgs = 100;

%% Test 1 - size of configurations

for h = 1:numel(HilbertArray)
    disp(['Testing Hilbert: ' HilbertIDs{h}]);
    Cfg_len_cond = TestCfgLength(HilbertArray{h});
    assert(Cfg_len_cond,'Test failed: configuration length');
end

%% Test 2 - difference self-consistency

for h = 1:numel(HilbertArray)
    disp(['Testing Hilbert: ' HilbertIDs{h}]);
    Diff_cond = TestPropMove(HilbertArray{h});
    assert(Diff_cond,'Test failed: difference consistency');
end

%% Test 3 - configuration restrictions

for h = 1:numel(HilbertArray)
    disp(['Testing Hilbert: ' HilbertIDs{h}]);
    [lower_cond, upper_cond, local_cond, bad_cfgs] = TestCfgBounds(HilbertArray{h},...
        N_cfgs,Cfg_lower_bounds(h),Cfg_upper_bounds(h));
    if ~(lower_cond && upper_cond && local_cond)
        disp(['Example bad configuration: ' num2str(bad_cfgs(1,:))]);
    end
    assert(lower_cond&&upper_cond&&local_cond,'Test failed: configuration restrictions')
end