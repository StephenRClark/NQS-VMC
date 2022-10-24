N = 6; NhV = 1:2*N;

FidVec = zeros(numel(NhV),1); EneVec = zeros(numel(NhV),1);
Fid0Vec = zeros(numel(NhV),1); Ene0Vec = zeros(numel(NhV),1);

AnsStr0 = 'Plus-NQSOH Nh'; AnsStrEnd = []; % ' FT'; % 

% SysStr = 'AKLT-SZ 1D'; load([SysStr ' N ' num2str(N) ' exact ground state.mat']);
% GS = gs_sz; H = H_sz; basis = basis_sz; EN = en_sz; 
SysStr = 'AKLT-XYZ 1D'; load([SysStr ' N ' num2str(N) ' exact ground state.mat']);
GS = gs_xyz; H = H_xyz; basis = basis_xyz; EN = en_xyz;

N0Inds = (abs(GS)>1e-15);

for n = 1:numel(NhV)
    AnsStr = [AnsStr0 ' ' num2str(NhV(n)) AnsStrEnd];
    load([SysStr ' N ' num2str(N) ' ' AnsStr ' Logs.mat']);
%     Psi = AnsatzObj.Modifier{1}.PsiGenerate(basis); Psi = Psi./sign(Psi(1)); 
%     Psi0 = Psi(N0Inds)/sqrt(sum(abs(Psi(N0Inds)).^2));
%     GS0 = GS(N0Inds)/sqrt(sum(abs(GS(N0Inds)).^2));
%     FidPsi = abs(GS'*Psi)^2; FidPsi0 = abs(GS0'*Psi0)^2; InfidPsi = 1 - FidPsi; InfidPsi0 = 1 - FidPsi0;
%     EnePsi = real(Psi'*H*Psi); EnePsi0 = real(Psi0'*H(N0Inds,N0Inds)*Psi0);
    FidVec(n) = FidPsi; Fid0Vec(n) = FidPsi0;    
    EneVec(n) = EnePsi; Ene0Vec(n) = EnePsi0;
%     save([SysStr ' N ' num2str(N) ' ' AnsStr ' Logs.mat'],...
%         'AnsatzObj','EnIter','RunTime','EnePsi','EnePsi0','FidPsi','FidPsi0','InfidPsi','InfidPsi0');
end

title('Energy versus hidden unit number');
figure(1); plot(NhV,EneVec,'-x'); hold on; plot(NhV,Ene0Vec,'-o');

figure(2); semilogy(NhV,1-FidVec,'-x'); hold on; semilogy(NhV,1-Fid0Vec,'-o');
title('Infidelity versus hidden unit number');
 
save([SysStr ' N ' num2str(N) ' ' AnsStr0 AnsStrEnd ' scan expectation values.mat'],'NhV',...
    'EneVec','FidVec','Ene0Vec','Fid0Vec');