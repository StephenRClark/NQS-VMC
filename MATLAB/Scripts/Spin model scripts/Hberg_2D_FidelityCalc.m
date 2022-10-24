L = 4; N = L^2; AlphaVec = [1 2 3 4 6 8 12 16];
load(['Heisenberg 2D N ' num2str(N) ' exact ground state.mat']);
FidA = zeros(numel(AlphaVec),1); EnePsiA = zeros(numel(AlphaVec),1);
EneErrorA = zeros(numel(AlphaVec),1);
for a = 1:numel(AlphaVec)
    AnsStr = ['Plus-NQSTI Alpha ' num2str(AlphaVec(a)) ' FT3'];
    load(['Heisenberg 2D N ' num2str(N) ' ' AnsStr ' Logs.mat'],'AnsatzObj');
    PsiA = PsiGenerate(AnsatzObj.Modifier{1},basis);
    FidA(a) = abs(PsiA'*gs)^2; EnePsiA(a) = PsiA'*H*PsiA;
    EneErrorA(a) = abs((EnePsiA(a) - gs_en)/gs_en);
end
figure(1); plot(AlphaVec,FidA,'-o'); hold on;
figure(2); plot(AlphaVec,EnePsiA,'-o'); hold on; plot(AlphaVec,gs_en*ones(numel(AlphaVec),1));
figure(3); plot(AlphaVec,EneErrorA,'-x'); hold on;