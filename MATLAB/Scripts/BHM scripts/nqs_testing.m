addpath(genpath("..\..\Ansatze"),genpath("..\..\Graphs"),genpath("..\..\Hilbert"));

Dim = [4 4]; N = prod(Dim); Bound = [1 1]; LVecs = eye(2);
GraphObj = HypCub(Dim,Bound,LVecs,1);

Nmax = 4; HilbertObj = Bose(N,N,Nmax);

ModParams.Alpha = 1; ModParams.HDim = 5; ModParams.AlphaP = 1;
ModParams.a = 0; ModParams.b = 0; ModParams.W = 0;
ModParams.nmag = 0; ModParams.nphs = 0;

NQS_test = NQS(HilbertObj,GraphObj,ModParams,1);
NQSA_test = NQSA(HilbertObj,GraphObj,ModParams,1);
NQSB_test = NQSB(HilbertObj,GraphObj,ModParams,1);
NQSU_test = NQSU(HilbertObj,GraphObj,ModParams,1);
NQSM_test = NQSM(HilbertObj,GraphObj,ModParams,1);

NewNQSArray = {NQS_test; NQSA_test; NQSB_test; NQSU_test; NQSM_test};

Cfg_test = HilbertObj.RandomCfg(); [Diff,CfgD_test] = HilbertObj.PropMove(Cfg_test);
RatioTester = [Cfg_test.occ.'; CfgD_test.occ.'];
CfgP_test = Cfg_test; Ncfgs = 100;
TestCfgMat = zeros(Ncfgs,N); TestCfgMat(1,:) = CfgP_test.occ.';
for n = 2:Ncfgs
    [~,CfgP_test] = HilbertObj.PropMove(CfgP_test);
    TestCfgMat(n,:) = CfgP_test.occ.';
end

for a = 1:numel(NewNQSArray)
TestPass = true;
% Testing ParamLoad and ParamList:
Params = linspace(0,1,NewNQSArray{a}.Np).';
NewNQSArray{a} = NewNQSArray{a}.ParamLoad(Params(:));
Params_read = NewNQSArray{a}.ParamList();
dParams = Params - Params_read;
if sum(abs(dParams))>1e-10
    disp('Test failed: ParamList');
    disp(['Original parameters: ' num2str(Params.')]);
    disp(['Listed parameters: ' num2str(Params_read.')]);
    TestPass = false;
end

% Testing PsiUpdate:
NewNQSArray{a} = NewNQSArray{a}.PsiUpdate(-Params);
Params_read = NewNQSArray{a}.ParamList();
if (sum(abs(Params_read)))>1e-15
    disp('Test failed: PsiUpdate');
    disp(['Resulting parameters: ' num2str(Params_read).']);
    TestPass = false;
end

NewNQSArray{a} = NewNQSArray{a}.ParamLoad(Params);

% Test that PrepPsi works correctly:
NewNQSArray{a} = NewNQSArray{a}.PrepPsi(Cfg_test);

% Testing PsiRatio works:
[Ratio_df, Update] = NewNQSArray{a}.PsiRatio(Diff);

% Testing PsiGenerate on small set:
Psi_gen = NewNQSArray{a}.PsiGenerate(RatioTester);
Ratio_wf = Psi_gen(2)/Psi_gen(1);
dRatio = Ratio_df - Ratio_wf;
if abs(dRatio/Ratio_wf) > 1e-10
    disp('Test failed: Ratio');
    disp(['Ratio from PsiGenerate: ' num2str(Ratio_wf)]);
    disp(['Ratio from PsiRatio: ' num2str(Ratio_df)]);
    TestPass = false;
end

% Testing dLogp gives correct number of values:
dLogp = NewNQSArray{a}.LogDeriv(Cfg_test);
if length(dLogp) ~= NewNQSArray{a}.Np
    disp('Test failed: dLogp');
    disp(['Length of dLogp: ' num2str(length(dLogp))]);
    disp(['Number of parameters: ' num2str(NewNQSArray{a}.Np)]);
    TestPass = false;
end

% Testing AddHidden works correctly:
NewNQSA2 = NewNQSArray{a}.AddHidden(ModParams);
Params_read = NewNQSA2.ParamList();
if (abs(sum(Params_read)-sum(Params))>1e-12) || (numel(Params_read)<numel(Params))
    disp('Test failed: AddHidden');
    disp(['Length of new Params: ' num2str(length(Params_read))]);
    disp(['New number of parameters: ' num2str(NewNQSA2.Np)]);
    disp(['Alpha of new NQS: ' num2str(NewNQSA2.Alpha)]);
    disp(['Size of OptInds: ' num2str(size(NewNQSA2.OptInds,1)) ' '...
        num2str(size(NewNQSA2.OptInds,2))]);
    TestPass = false;
end

if TestPass
    disp(['All tests passed for ' class(NewNQSArray{a})]);
else 
    disp(['Some tests failed for ' class(NewNQSArray{a})]);
end

end