% Script for constructing and optimising Ansatz objects with number-hidden
% NQS Modifiers and Bose condensate References for the Bose Hubbard model
% on a 1D lattice with periodic boundary conditions.

L = 8; N = L^2; Nmax = 4; % HilbertObj = Bose(N,N,Nmax); GraphObj = HypCub([L L],[1 1],eye(2),1);

NStr = ' N '; U = [1 4 8 12 14 15 16 17 18 19 20 21 22 24 28 32 40 48];

ModParams.Nh = N; ModParams.Alpha = 1; ModParams.nmag = 0.05; ModParams.nphs = 0; 
ModParams.a = -0.01; ModParams.b = 0.01; ModParams.W = 0.01; ModParams.X = -0.01;
ModParams.A = -0.1; ModParams.B = -0.01;

% AnsStr0 = 'BECR-Gutz';
% AnsStr0 = ['BECR-NQSSHTI-HD' num2str(Nmax+1) '-bB Alpha 1']; % Ansatz name string for easy identification later.
AnsStr =  ['BECR-NQSSHTI-HD2-AW Alpha 1'];

urange = 1:numel(U);

for u = urange
    load(['BHM 2D U ' num2str(U(u)) NStr num2str(N) ' ' AnsStr ' Logs.mat'],'AnsatzObj');
    % NQSObj = AnsatzObj.Modifier{1}; NQSObj = NQSObj.AddHidden(ModParams);
    % NQSObj.OptInds(1) = 0; AnsatzObj = AnsatzObj.ModReplace(NQSObj,1);
    % NQSObj = NQSMHTI(HilbertObj,GraphObj,ModParams,1); NQSObj.OptInds(1) = 0; 
    % dP = zeros(NQSObj.Np,1); dP(2) = -AnsatzObj.Modifier{1}.G; 
    % dP(5) = -0.1; dP(5+N) = 0.1; NQSObj = PsiUpdate(NQSObj,dP);
    % AnsatzObj = ModReplace(AnsatzObj,NQSObj,1); 
    AnsProp = AnsatzObj.PropertyList;
    AnsProp.Modifier{1}.Graph.Bound = [1 1]; AnsProp.Reference.Graph.Bound = [1 1];
    save(['BHM 2D U ' num2str(U(u)) NStr num2str(N) ' ' AnsStr ' setup.mat'],'AnsProp');
end

txtstr = 'nqsbh';

fileID = fopen(['bhm_' num2str(N) '_' txtstr '_setup.txt'],'w');
for u = urange
    fprintf(fileID,['BHM 2D U ' num2str(U(u)) NStr num2str(N) ' ' AnsStr ' setup.mat']);
    if u < urange(end)
        fprintf(fileID,'\r\n');
    end
end
fclose(fileID);