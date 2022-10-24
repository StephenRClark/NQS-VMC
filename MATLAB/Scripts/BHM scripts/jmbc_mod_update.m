UVec = [1 16 23 32]; N = 100; NStr = [' N ' num2str(N)];

Nmax = 4; HilbertObj = Bose(N,N,Nmax);

AnsStrArray = {'BECR-JMBC-NQSA'; 'BECR-JMBC-NQSB-HDim5'; 'BECR-JMBC-NQSU-VDim5'};

for a = 1:numel(AnsStrArray)
    for u = 1:numel(UVec)
        U = UVec(u);
        load(['BHM 2D U ' num2str(U) NStr ' ' AnsStrArray{a} ' Logs.mat']);
        N_mod = numel(AnsatzObj.Modifier);
        for n = 1:N_mod
            Mod = AnsatzObj.Modifier{n};
            if size(Mod.OptInds,2)==1
                Mod.OptInds = [Mod.OptInds, Mod.OptInds];
                AnsatzObj = AnsatzObj.ModReplace(Mod,n);
            end
        end
        save(['BHM 2D U ' num2str(U) NStr ' ' AnsStrArray{a} ' Logs.mat'],'AnsatzObj',...
                'BiBj','DbHl','DiDj','EneGS','EnIter','EvalTime','HiHj','NiNj',...
                'OcFr','Params','RunDate','RunTime','VarE','VarN');
    end
end