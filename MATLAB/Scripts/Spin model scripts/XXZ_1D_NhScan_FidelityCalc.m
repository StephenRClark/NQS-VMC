DV = [0 1]; Nvec{1} = 1:20; Nvec{2} = 1:20;

for d = 1:numel(DV)
    
    Delta = DV(d); N = 20;
    
    NhV = Nvec{d};
    
    FidVec = zeros(numel(NhV),1); EneVec = zeros(numel(NhV),1); EneErrVec = zeros(numel(NhV),1);
    Fid0Vec = zeros(numel(NhV),1); Ene0Vec = zeros(numel(NhV),1); Ene0ErrVec = zeros(numel(NhV),1);
    
    AnsStr0 = 'Plus-NQS'; AnsStrEnd = ' JSW';
    
    load(['XXZ Delta ' num2str(Delta) ' N ' num2str(N) ' exact ground state.mat']);
    GS = gs_ms; HMT = H_ms; Sz0Inds = (sum(basis,2)==0);
    
    for n = 1:numel(NhV)
        AnsStr = [AnsStr0 ' Nh ' num2str(NhV(n)) AnsStrEnd];
        load(['XXZ 1D Jz ' num2str(Delta) ' N ' num2str(N) ' ' AnsStr ' Logs.mat']);
        Psi = AnsatzObj.Modifier{1}.PsiGenerate(basis); Psi = Psi./sign(Psi(1));
        Psi0 = Psi(Sz0Inds)/sqrt(sum(abs(Psi(Sz0Inds)).^2));
        GS0 = GS(Sz0Inds)/sqrt(sum(abs(GS(Sz0Inds)).^2));
        FidPsi = abs(GS'*Psi)^2; FidVec(n) = FidPsi; FidPsi0 = abs(GS0'*Psi0)^2; Fid0Vec(n) = FidPsi0;
        EnePsi = real(Psi'*HMT*Psi); EneVec(n) = EnePsi; EnePsi0 = real(Psi0'*HMT(Sz0Inds,Sz0Inds)*Psi0); Ene0Vec(n) = EnePsi0;
        EneErr = abs((EnePsi - en)/en); EneErrVec(n) = EneErr; Ene0Err = abs((EnePsi0 - en)/en); Ene0ErrVec(n) = Ene0Err;
        %     save(['XXZ 1D Jz ' num2str(Delta) ' N ' num2str(N) ' ' AnsStr ' Logs.mat'],...
        %         'AnsatzObj','EnIter','RunTime','FidPsi','FidPsi0',...
        %         'EnePsi','EnePsi0','EneErr','Ene0Err'); % ,'Psi','Psi0');
    end
    
    figure(1); plot(NhV,EneVec,'-x'); hold on; plot(NhV,Ene0Vec,'-o');
    title('Energy versus hidden unit number');
    
    figure(2); plot(NhV,FidVec,'-x'); hold on; plot(NhV,Fid0Vec,'-o');
    title('Fidelity versus hidden unit number');
    
    figure(3); plot(NhV,EneErrVec,'-x'); hold on; plot(NhV,Ene0ErrVec,'-o');
    title('Energy error versus hidden unit number');
    
    save(['XXZ 1D Jz ' num2str(Delta) ' N ' num2str(N) ' ' AnsStr0 AnsStrEnd ' Nh scan expectation values.mat'],'NhV',...
        'EneVec','FidVec','EneErrVec','Ene0Vec','Fid0Vec','Ene0ErrVec');
end

% % Reshuffling lines for sorting the raw values.
w_p = w_d0; w_p0 = w_p;
[Wm,Im] = max(abs(w_p),[],2); [Ws,Is] = sort(Wm,'descend');
for h = 1:20
    Wsign = sign(w_p(Is(h),Im(Is(h)))); w_p0(h,:) = Wsign*circshift(w_p(Is(h),:),-Im(Is(h))+h);
end
w_d0_rs = w_p0;
