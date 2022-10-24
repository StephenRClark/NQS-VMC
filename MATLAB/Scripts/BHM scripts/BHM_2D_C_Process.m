clear; clc;

L = 6; N = L^2; Dim = [L L];

U = [1 4 8 12 14 15 16 17 18 19 20 21 22 24 28 32 40 48];

% AnsStr = 'BECR-NQSNHTI-HD5-bB Alpha 1';
% AnsStr = 'BECR-NQSNHTI-HD5-NNMB Alpha 2';
% AnsStr = 'BECR-NQSSHTI-HD5-bB Alpha 1';
% AnsStr = 'BECR-NQSSHTI-HD5 Alpha 2';
% AnsStr = 'BECR-Jast'; 
% AnsStr = 'BECR-Jast-NNMB';
% AnsStr = 'BECR-NQSSHTI-HD2-AW Alpha 1';
% AnsStr = 'BECR-NQSSHTI-HD2-AW Alpha 2';
% AnsStr = 'BECR-NQSMHTI-HD5-LR Alpha 1';
% AnsStr = 'BECR-NQSOHTI-VD5 Alpha 1';

urange = 1:numel(U); % (1:3) + 3; %  % Run this first
s = 5;
for u = urange
    load(['BHM 2D U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' setup.mat']); %   ' num2str(s) '
    save(['BHM 2D U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' setup ' num2str(s) '.mat'],'AnsProp','StartParams',...
        'ParamPass1','EnIter1','ParamPass2','EnIter2','ParamPass3','EnIter3',...
        'RunTime','EvalTime','EneGS','EnSq','VarN','NiNj','DbHl','OcFr','BiBj');
    delete(['BHM 2D U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' setup.mat']);
    EnIter = [EnIter1; EnIter2; EnIter3]; EneTotal = EneGS; EneGS = EneGS / N;
    ParamSets = [StartParams ParamPass1 ParamPass2 ParamPass3];
    save(['BHM 2D U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' C results ' num2str(s) '.mat'],...
        'EnIter','AnsProp','RunTime','ParamSets','EvalTime',...
        'EneGS','EneTotal','EnSq','VarN','NiNj','DbHl','OcFr','BiBj');
    figure(u); plot(real(EnIter));
end