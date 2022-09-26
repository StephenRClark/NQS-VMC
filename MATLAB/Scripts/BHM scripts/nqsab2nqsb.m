UVec = [1 16 23 32]; N = 100; NStr = [' N ' num2str(N)];

Nmax = 4; HilbertObj = Bose(N,N,Nmax);

ModParams.Alpha = 2; ModParams.HDim = 5;
ModParams.a = 0; ModParams.b = 0; ModParams.W = 0;
ModParams.nmag = 0; ModParams.nphs = 0;

AnsStrOld = {'BECR-NQS-AB-HDim5 Alpha 2'};

AnsStrNew = {'BECR-NQSB-HDim5 Alpha 2'};

% Testing setup
TestPass = true; Ncfgs = 100;
TestCfg = HilbertObj.RandomCfg(); TestCfgP = TestCfg; [Diff,TestCfgD] = HilbertObj.PropMove(TestCfg);
TestCfgMat = zeros(Ncfgs,N); TestCfgMat(1,:) = TestCfg.occ.';
for n = 2:Ncfgs
    [~,TestCfgP] = HilbertObj.PropMove(TestCfgP);
    TestCfgMat(n,:) = TestCfgP.occ.';
end

for a = 1:numel(AnsStrOld)
    for u = 1:numel(UVec)
        TestPass = true;
        U = UVec(u);
        load(['BHM 2D U ' num2str(U) NStr ' ' AnsStrOld{a} ' Logs.mat']);
        NQSObjOld = AnsatzObj.Modifier{1}; ParamsOld = NQSObjOld.ParamList;
        % Setup new Modifier.
        Params = ParamsOld; ParamCap = 10;
        % Perform conversion for spin-hidden versions.
        if isa(NQSObjOld,'NQSSHTI')
            % Require rescaling of hidden bias to accommodate linear shift
            b_old = ParamsOld([3 4]); B_old = ParamsOld([5 6]); 
            b_new = b_old - 4*B_old; Params([3 4]) = b_new;
            while max(abs(b_new)) > ParamCap
                ParamCap = ParamCap + 1;
            end
        end
        OptIndsOld = NQSObjOld.OptInds; OptIndsNew = [OptIndsOld, (OptIndsOld*0)];
        NQSObjNew = NQSB(HilbertObj, NQSObjOld.Graph, ModParams, 1); NQSObjNew.ParamCap = ParamCap;
        NQSObjNew.OptInds = OptIndsNew; NQSObjNew = NQSObjNew.ParamLoad(Params);
        
        % Testing basic quantities to ensure equivalence
        NQSObjOld = NQSObjOld.PrepPsi(TestCfg); NQSObjNew = NQSObjNew.PrepPsi(TestCfg);
        % Theta check:
        ThetaOld = NQSObjOld.Theta; ThetaNew = NQSObjNew.Theta;
        if isa(NQSObjOld,'NQSSHTI')
            % Account for b shift
            ThetaNew(1:N) = ThetaNew(1:N) + 4*B_old(1);
            ThetaNew((1:N)+N) = ThetaNew((1:N)+N) + 4*B_old(2);
        end        
        dTheta = (ThetaNew - ThetaOld)./ThetaOld;
        if sum(abs(dTheta))>1e-10
            disp('Theta check failed.'); TestPass = false;
            disp(['Summed Theta difference: ' num2str(sum(abs(dTheta)))]);
        end
        % Ratio check:
        RatioOld = NQSObjOld.PsiRatio(Diff); RatioNew = NQSObjNew.PsiRatio(Diff);
        dRatio = (RatioNew - RatioOld)/RatioOld;
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
            AnsatzObj = AnsatzObj.ModReplace(NQSObjNew,1);
            save(['BHM 2D U ' num2str(U) NStr ' ' AnsStrNew{a} ' Logs.mat'],'AnsatzObj',...
                'BiBj','DbHl','DiDj','EneGS','EnIter','EvalTime','HiHj','NiNj',...
                'OcFr','Params','RunDate','RunTime','VarE','VarN');
        else
            disp(['Test failure at U = ' num2str(U)]);
        end
    end
end
