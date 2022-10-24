Nt = 4; N = 3*Nt; NhV = 8;

FidVec = zeros(numel(NhV),1); Fid0Vec = zeros(numel(NhV),1);

AnsStr0 = 'Plus-NQS Nh'; AnsStrEnd = []; % ' FT'; % 

load(['AKLT-XYZ 1D N ' num2str(Nt) ' exact ground state.mat']);

GS = gs_xyz; basis = basis_u; N0Inds = (abs(GS)>1e-15); 

for n = 1:numel(NhV)
    AnsStr = [AnsStr0 ' ' num2str(NhV(n)) AnsStrEnd];
    load(['GS 1D N ' num2str(N) ' ' AnsStr ' Logs.mat']);
    Psi = AnsatzObj.Modifier{1}.PsiGenerate(basis); Psi = Psi/sign(Psi(1)); 
    Psi0 = Psi(N0Inds)/sqrt(sum(abs(Psi(N0Inds)).^2));
    GS0 = GS(N0Inds)/sqrt(sum(abs(GS(N0Inds)).^2));
    FidPsi = abs(GS'*Psi)^2; FidVec(n) = FidPsi; FidPsi0 = abs(GS0'*Psi0)^2; Fid0Vec(n) = FidPsi0;
    save(['GS 1D N ' num2str(N) ' ' AnsStr ' Logs.mat'],...
        'AnsatzObj','EnIter','RunTime','FidPsi','FidPsi0');
end

figure(1); plot(NhV,FidVec,'-x'); hold on; plot(NhV,Fid0Vec,'-o');
title('Unary subspace fidelity versus hidden unit number');

% save(['GS 1D N ' num2str(N) ' ' AnsStr0 AnsStrEnd ' scan expectation values.mat'],'NhV',...
%     'FidVec','Fid0Vec');