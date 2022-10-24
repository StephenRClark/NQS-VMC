UVec = [1 16 23 32]; N = 100; NStr = [' N ' num2str(N)];

Nmax = 4; HilbertObj = Bose(N,N,Nmax);

ModParams.Alpha = 1; ModParams.VDim = 5;
ModParams.a = 0; ModParams.b = 0; ModParams.W = 0;
ModParams.nmag = 0; ModParams.nphs = 0;

AnsStrOld = {'BECR-JHD-NQS-U-VDim5'};

AnsStrNew = {'BECR-JMBC-NQSU-VDim5'};

% Testing setup
TestPass = true; Ncfgs = 100;
TestCfg = HilbertObj.RandomCfg(); TestCfgP = TestCfg; [Diff,TestCfgD] = HilbertObj.PropMove(TestCfg);
TestCfgMat = zeros(Ncfgs,N); TestCfgMat(1,:) = TestCfg.occ.';
for n = 2:Ncfgs
    [~,TestCfgP] = HilbertObj.PropMove(TestCfgP);
    TestCfgMat(n,:) = TestCfgP.occ.';
end

for a = 1%:numel(AnsStrOld)
    for u = 1:numel(UVec)
        U = UVec(u);
        load(['BHM 2D U ' num2str(U) NStr ' ' AnsStrOld{a} ' Logs.mat']);
        NQSObjOld = AnsatzObj.Modifier{3}; ParamsOld = NQSObjOld.ParamList;
        % Separate out parameters corresponding to v = 0:
        a_old = ParamsOld(1:(Nmax+1)); b_old = ParamsOld(Nmax+2);
        w_old = reshape(ParamsOld((1:((Nmax+1)*N))+Nmax+2),Nmax+1,N);
        a_0 = a_old(1); w_0 = w_old(1,:);
        a_new = ParamsOld((1:Nmax)+1) - a_0;
        b_new = b_old + sum(w_0);
        w_new = reshape((w_old(2:end,:)-w_0),Nmax*N,1);
        % Setup new Modifier.
        Params = [a_new; b_new; w_new];
        OptIndsOld = NQSObjOld.OptInds; OptIndsNew = [OptIndsOld, (OptIndsOld*0)];
        NQSObjNew = NQSU(HilbertObj, NQSObjOld.Graph, ModParams, 1); NQSObjNew.ParamCap = 10;
        NQSObjNew.OptInds = OptIndsNew; NQSObjNew = NQSObjNew.ParamLoad(Params);
        
        % Testing basic quantities to ensure equivalence
        NQSObjOld = NQSObjOld.PrepPsi(TestCfg); NQSObjNew = NQSObjNew.PrepPsi(TestCfg);
        % Theta check:
        ThetaOld = NQSObjOld.Theta; ThetaNew = NQSObjNew.Theta;
        dTheta = ThetaNew - ThetaOld;
        if sum(abs(dTheta))>1e-10
            disp('Theta check failed.'); TestPass = false;
            disp(['Summed Theta difference: ' num2str(sum(abs(dTheta)))]);
        end
        % Ratio check:
        RatioOld = NQSObjOld.PsiRatio(Diff); RatioNew = NQSObjNew.PsiRatio(Diff);
        dRatio = abs((RatioNew - RatioOld)/RatioOld);
        if abs(dRatio)>1e-10
            disp('Ratio check failed.'); TestPass = false;
            disp(['Ratio difference: ' num2str(RatioNew-RatioOld)]);
            disp(['Ratio old: ' num2str(RatioOld) ' vs Ratio new: ' num2str(RatioNew)]);
        end
        % Psi check:
        PsiOld = NQSObjOld.PsiGenerate(TestCfgMat); PsiNew = NQSObjNew.PsiGenerate(TestCfgMat);
        Infid = 1 - abs((PsiOld.'*PsiNew)^2);
        if Infid > 1e-10
            disp('Infidelity check failed.'); TestPass = false;
            disp(['Fidelity: ' num2str(1 - Infid)]);
            disp(['Infidelity: ' num2str(Infid)]);
        end
        
        if TestPass
            disp(['Success for U = ' num2str(U)]);
            AnsatzObj = AnsatzObj.ModReplace(NQSObjNew,3);
            save(['BHM 2D U ' num2str(U) NStr ' ' AnsStrNew{a} ' Logs.mat'],'AnsatzObj',...
                'BiBj','DbHl','DiDj','EneGS','EnIter','EvalTime','HiHj','NiNj',...
                'OcFr','Params','RunDate','RunTime','VarE','VarN');
            else
            disp(['Test failure at U = ' num2str(U)]);
        end
    end
end