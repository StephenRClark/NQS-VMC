addpath(genpath("../Ansatze"),genpath("../Graphs"),genpath("../Hilbert"));

% Testing with 4x4 grid of bosons.

L = 4; Dim = [L L]; N = prod(Dim); Bound = [1 1]; LVecs = eye(2); 
GraphObj = HypCub(Dim,Bound,LVecs,1); 

Nmax = 4; HilbertObj = Bose(N,N,Nmax);

ModParams.Alpha = 1; ModParams.HDim = 5; ModParams.AlphaP = 1;
ModParams.a = 0.2; ModParams.b = -0.5; ModParams.W = 0.1;
ModParams.nmag = 0.01; ModParams.nphs = 0;

Cfg_test = HilbertObj.RandomCfg(); [Diff,CfgD_test] = HilbertObj.PropMove(Cfg_test);
RatioTester = [HilbertObj.FullCfg(Cfg_test).'; HilbertObj.FullCfg(CfgD_test).'];
CfgP_test = Cfg_test; Ncfgs = 100;
TestCfgMat = zeros(Ncfgs,N); TestCfgMat(1,:) = HilbertObj.FullCfg(CfgP_test).';
for n = 2:Ncfgs
    [~,CfgP_test] = HilbertObj.PropMove(CfgP_test);
    TestCfgMat(n,:) = HilbertObj.FullCfg(CfgP_test);
end

% Generate set of NQS with lossless conversions:
NQSAObj = NQSA(HilbertObj,GraphObj,ModParams,1);
NQSBObj = NQSB(HilbertObj,GraphObj,ModParams,1);
NQSMObj = NQSM(HilbertObj,GraphObj,ModParams,1);

% New set of converted NQS from lossless conversions:
NQSA2U = nqsa2nqsu(NQSAObj,HilbertObj);
NQSB2C = nqsb2nqsc(NQSBObj,HilbertObj);
NQSM2U = nqsm2nqsu(NQSMObj,HilbertObj);

OldNQSArray = {NQSAObj; NQSBObj; NQSMObj};
NewNQSArray = {NQSA2U; NQSB2C; NQSM2U};

%% Test 1: PsiRatio
for a = 1:numel(OldNQSArray)
    disp(['Old Modifier class: ' class(OldNQSArray{a})]);
    disp(['New Modifier class: ' class(NewNQSArray{a})]);
    OldNQSArray{a} = OldNQSArray{a}.PrepPsi(Cfg_test);
    NewNQSArray{a} = NewNQSArray{a}.PrepPsi(Cfg_test);
    PsiRatioOld = OldNQSArray{a}.PsiRatio(Diff);
    PsiRatioNew = NewNQSArray{a}.PsiRatio(Diff);
    Psi_rat_cond = abs(1-(PsiRatioNew/PsiRatioOld))<1e-10;
    assert(Psi_rat_cond,'Test failed: PsiRatio')
    if ~Psi_rat_cond
        disp(['Test failure for conversion between ' class(OldNQSArray{a}) ...
            ' to ' class(NewNQSArray{a})]);
    end
end
%% Test 2: PsiGenerate (medium set)
for a = 1:numel(OldNQSArray)
    disp(['Old Modifier class: ' class(OldNQSArray{a})]);
    disp(['New Modifier class: ' class(NewNQSArray{a})]);
    PsiOld = OldNQSArray{a}.PsiGenerate(TestCfgMat);
    PsiNew = NewNQSArray{a}.PsiGenerate(TestCfgMat);
    Fid = PsiNew'*PsiOld; Psi_fid_cond = (1-Fid)<1e-10;
    assert(Psi_fid_cond,'Test failed: PsiGenerate');
    if ~Psi_fid_cond
        disp(['Test failure for conversion between ' class(OldNQSArray{a}) ...
            ' to ' class(NewNQSArray{a})]);
        disp(['Infidelity: ' num2str(1-Fid)']);
    end
end