addpath(genpath("..\Ansatze"),genpath("..\Graphs"),genpath("..\Hilbert"));
addpath("Test functions/Modifiers");

% Testing with 4x4 grid of spin-1s.

Dim = [4 4]; N = prod(Dim); Bound = [1 1]; LVecs = eye(2);
GraphObj = HypCub(Dim,Bound,LVecs,1);

S = 1; SzT = 0; HilbertObj = Spin(N,S,SzT);

ModParams.Alpha = 1; ModParams.HDim = 5; ModParams.AlphaP = 1;
ModParams.a = 0; ModParams.b = 0; ModParams.W = 0;
ModParams.nmag = 0; ModParams.nphs = 0; ModParams.Js = 0;

% Create Modifiers to test their core functions:

NQS_test = NQS(HilbertObj,GraphObj,ModParams,1);
NQSU_test = NQSU(HilbertObj,GraphObj,ModParams,1);
NQSS1_test = NQSS1(HilbertObj,GraphObj,ModParams,1);
Jast_test = Jast(HilbertObj,GraphObj,ModParams,1); Jast_test.NormFlag = 0;

ModifierArray = {NQS_test; NQSU_test; NQSS1_test; Jast_test};

Cfg_test = HilbertObj.RandomCfg(); [Diff,CfgD_test] = HilbertObj.PropMove(Cfg_test);
RatioTester = [HilbertObj.FullCfg(Cfg_test).'; HilbertObj.FullCfg(CfgD_test).'];
CfgP_test = Cfg_test; Ncfgs = 100;
TestCfgMat = zeros(Ncfgs,N); TestCfgMat(1,:) = HilbertObj.FullCfg(CfgP_test).';
for n = 2:Ncfgs
    [~,CfgP_test] = HilbertObj.PropMove(CfgP_test);
    TestCfgMat(n,:) = HilbertObj.FullCfg(CfgP_test);
end

ParamsArray = cell(numel(ModifierArray),1);
for a = 1:numel(ModifierArray)
    ParamsArray{a} = linspace(0,1,ModifierArray{a}.Np).';
end

%% Test 1: ParamLoad and ParamList
for a = 1:numel(ModifierArray)
    disp(['Testing Modifier class: ' class(ModifierArray{a})]);
    Param_list_cond = TestParams_Load_List(ModifierArray{a},ParamsArray{a});
    assert(Param_list_cond,'Test failed: ParamLoad / ParamList');
    if ~Param_list_cond
        disp(['Test failure for class: ' class(ModifierArray{a})]);
    end
end
%% Test 2: PsiUpdate
for a = 1:numel(ModifierArray)
    disp(['Testing Modifier class: ' class(ModifierArray{a})]);
    Psi_update_cond = TestPsiUpdate(ModifierArray{a});
    assert(Psi_update_cond,'Test failed: PsiUpdate');
    if ~Psi_update_cond
        disp(['Test failure for class: ' class(ModifierArray{a})]);
    end
end
%% Test 3: PsiRatio and PsiGenerate (small set)
for a = 1:numel(ModifierArray)
    disp(['Testing Modifier class: ' class(ModifierArray{a})]);
    Psi_rat_cond = TestPsiRatio(ModifierArray{a},Cfg_test,Diff,RatioTester);
    assert(Psi_rat_cond,'Test failed: PsiRatio / PsiGenerate')
    if ~Psi_rat_cond
        disp(['Test failure for class: ' class(ModifierArray{a})]);
    end
end
%% Test 4: PsiGenerate (medium set)
for a = 1:numel(ModifierArray)
    disp(['Testing Modifier class: ' class(ModifierArray{a})]);
    [Psi_len_cond,Psi_norm_cond] = TestPsiGenerate(ModifierArray{a},TestCfgMat);
    assert(Psi_len_cond,'Test failed: PsiGenerate (length)');
    assert(Psi_norm_cond,'Test failed: PsiGenerate (normalisation)');
    if ~Psi_len_cond || ~Psi_norm_cond
        disp(['Test failure for class: ' class(ModifierArray{a})]);
    end
end
%% Test 5: dLogp
for a = 1:numel(ModifierArray)
    disp(['Testing Modifier class: ' class(ModifierArray{a})]);
    ModifierArray{a} = ModifierArray{a}.PrepPsi(Cfg_test);
    dLogp_cond = TestdLogp(ModifierArray{a},Cfg_test);
    assert(dLogp_cond,'Test failed: dLogp');
    if ~dLogp_cond
        disp(['Test failure for class: ' class(ModifierArray{a})]);
    end
end
%% Test 6: AddHidden
for a = 1:numel(ModifierArray)
    disp(['Testing Modifier class: ' class(ModifierArray{a})]);
    if contains(class(ModifierArray{a}),'NQS') % Additional function to test with NQS
        [Params_sum_cond,Params_len_cond] = TestAddHidden(ModifierArray{a},ModParams);
        assert(Params_sum_cond,'Test failed: AddHidden (parameter sum)');
        assert(Params_len_cond,'Test failed: AddHidden (new parameter number)');
    else % Automatic pass for other Modifiers since this is not a relevant function.
        disp('AddHidden not a relevant method for this Modifier.')
        Params_sum_cond = true; Params_len_cond = true;
    end
    if ~Params_sum_cond || ~Params_len_cond
        disp(['Test failure for class: ' class(ModifierArray{a})]);
    end
end