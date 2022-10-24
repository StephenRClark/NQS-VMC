L = 6; N = L^2; Dim = [L L];

U = [1 4 8 12 14 15 16 17 18 19 20 21 22 24 28 32 40 48]; 

% Create objects to load parameters into.
HilbertObj = Bose(N,N,4); GraphObj = HypCub(Dim,[1 1],eye(2),1);
HilbertBH = Bose(N,N,1);

RefParams.SPH = - Graph2Array(GraphObj,0); Ref = BECR(HilbertObj,GraphObj,RefParams);

ModParams.nmag = 0; ModParams.nphs = 0;
ModParams.a = 0; ModParams.b = 0; ModParams.W = 0;

% AnsStr = 'BECR-NQSNHTI-HD5-bB Alpha 1'; ModParams.Alpha = 1; Mod = {NQSNHTI(HilbertObj,GraphObj,ModParams,1)};
% AnsStr = 'BECR-NQSNHTI-HD5-NNMB Alpha 2'; ModParams.Alpha = 1; ModParams.GMB = 0; Mod = {NQSNHTI(HilbertObj,GraphObj,ModParams,1); NNMB(HilbertObj,GraphObj,ModParams,1)};
AnsStr = 'BECR-NQSSHTI-HD5-bB Alpha 1'; ModParams.Alpha = 1; Mod = {NQSSHTI(HilbertObj,GraphObj,ModParams,1)};
% AnsStr = 'BECR-NQSSHTI-HD5 Alpha 2'; ModParams.Alpha = 2; Mod = {NQSSHTI(HilbertObj,GraphObj,ModParams,1)};
% AnsStr = 'BECR-Jast'; ModParams.Js = 0; Mod = {Jast(HilbertObj,GraphObj,ModParams,1)};
% AnsStr = 'BECR-Jast-NNMB'; ModParams.Js = 0; ModParams.GMB = 0; Mod = {Jast(HilbertObj,GraphObj,ModParams,1);NNMB(HilbertObj,GraphObj,ModParams,1)};
% AnsStr = 'BECR-NQSSHTI-HD2-AW Alpha 1'; ModParams.Alpha = 1; Mod = {NQSSHTI(HilbertBH,GraphObj,ModParams,1)};
% AnsStr = 'BECR-NQSSHTI-HD2-AW Alpha 2'; ModParams.Alpha = 1; Mod = {NQSSHTI(HilbertBH,GraphObj,ModParams,1)};
% AnsStr = 'BECR-NQSMHTI-HD5-LR Alpha 1'; ModParams.Alpha = 1; Mod = {NQSMHTI(HilbertObj,GraphObj,ModParams,1)};
% AnsStr = 'BECR-NQSOHTI-VD5 Alpha 1';

urange = 1:numel(U);

% Use this block to generate AnsatzObj from AnsProp.
for u = urange
    load(['BHM 2D U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' C results.mat']); % 
    AnsatzObj = Ansatz(Ref,Mod,HilbertObj); AnsatzObj = AnsatzObj.ParamLoad(ParamSets(:,4));
    VarE = (EnSq - (EneTotal)^2)/(N^2);
    save(['BHM 2D U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' C Logs.mat'],...
        'EnIter','AnsatzObj','RunTime','EvalTime','ParamSets',...
        'EneGS','EnSq','VarN','VarE','NiNj','DbHl','OcFr','BiBj','EvalTime','RunTime');
end

% % Use this block to generate AnsProp from AnsatzObj.
% for u = urange
%     load(['BHM 2D U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat'],'AnsatzObj'); % 
%     AnsProp = AnsatzObj.PropertyList; GraphProp = GraphObj.PropertyList;
%     AnsProp.Reference.Graph = GraphProp; AnsProp.Modifier{1}.Graph = GraphProp;
%     AnsProp.Modifier{1}.ParamCap = 10; % AnsProp.Modifier = {AnsProp.Modifier{1};NNMBProp};
%     save(['BHM 2D U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' setup.mat'],...
%         'AnsProp');
% end