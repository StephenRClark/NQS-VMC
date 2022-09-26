U = [1 4 8 12 14 15 16 17 18 19 20 21 22 24 28 32 40 48];

% AnsStr = 'BECR-NQSNHTI-HD5-bB Alpha 1'; N = 100; SetRange = 0:5; % N = 36; SetRange = 0:5; % N = 64; SetRange = 0:5; % 
% AnsStr = 'BECR-NQSNHTI-HD5-NNMB Alpha 2'; N = 36; SetRange = [0:3 5]; % N = 64; SetRange = 0:4; % N = 100; SetRange = 0:4; % 
% AnsStr = 'BECR-Jast'; N = 100; SetRange = 0:5; % N = 36; SetRange = 0:5; % N = 64; SetRange = 0:5; % 
% AnsStr = 'BECR-Jast-NNMB'; N = 100; SetRange = 0:5; % N = 36; SetRange = 0:5; % N = 64; SetRange = 0:5; % 
% AnsStr = 'BECR-NQSSHTI-HD2-AW Alpha 1'; N = 100; SetRange = 0:5; % N = 36; SetRange = 0:5; % N = 64; SetRange = 0:5; % 
% AnsStr = 'BECR-NQSSHTI-HD2-AW Alpha 2'; N = 100; SetRange = 0:5; % N = 36; SetRange = 0:5; % N = 64; SetRange = 0:5; % 
AnsStr = 'BECR-NQSSHTI-HD5-bB Alpha 1'; N = 100; SetRange = 0:5; % N = 36; SetRange = 0:5; % N = 64; SetRange = 0:5; % 
% AnsStr = 'BECR-NQSSHTI-HD5 Alpha 2'; N = 36; SetRange = 0:4; % N = 64; SetRange = 0:4; % N = 100; SetRange = 0:4; % 
% AnsStr = 'BECR-NQSMHTI-HD5-LR Alpha 1'; N = 64; SetRange = 0:3; % N = 100; SetRange = 0:3; % N = 36; SetRange = 0:3; % 

urange = 1:numel(U); 

EneList = zeros(numel(urange),numel(SetRange));

for s = 1:numel(SetRange)
    for u = urange
        load(['BHM 2D U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' C results ' num2str(SetRange(s)) '.mat'],'EneGS');
        EneList(u,s) = EneGS;
    end
end

RunDate = date;

[EneGSU,MinSets] = min(EneList,[],2); plot(U,EneGSU,'-o'); hold on;

save(['BHM 2D U Scan N ' num2str(N) ' ' AnsStr ' C best energy values.mat'],'EneGSU','MinSets','SetRange','urange','RunDate');

for u = 1:numel(U)
load(['BHM 2D U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' C results ' num2str(SetRange(MinSets(u))) '.mat']);
if numel(AnsProp.Modifier) == 2  
    AnsProp.Modifier{2}.Params = ParamSets(end,4); AnsProp.Modifier{1}.Params = ParamSets(1:(end-1),4);
elseif numel(AnsProp.Modifier) == 1
    AnsProp.Modifier{1}.Params = ParamSets(:,4);
end
save(['BHM 2D U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' setup.mat'],'AnsProp');
save(['BHM 2D U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' C results.mat'],'AnsProp','BiBj','DbHl','EnIter','EnSq',...
    'EneGS','EneTotal','EvalTime','NiNj','OcFr','ParamSets','RunTime','VarN');
end