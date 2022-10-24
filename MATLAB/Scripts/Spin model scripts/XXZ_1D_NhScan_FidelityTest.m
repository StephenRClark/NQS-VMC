Delta = 0; N = 20; load(['XXZ Delta ' num2str(Delta) ' N ' num2str(N) ' exact ground state.mat']);

AnsStr0 = 'Plus-NQS Nh '; AnsStrEnd = ' JSW FT'; Nh = 13; AnsStr = [AnsStr0 num2str(Nh) AnsStrEnd];
load(['XXZ 1D Jz ' num2str(Delta) ' N ' num2str(N) ' ' AnsStr ' Logs.mat'],'AnsatzObj');

NhV = 1:Nh;

FidExVec = zeros(numel(NhV),1); EneExVec = zeros(numel(NhV),1); EneErrExVec = zeros(numel(NhV),1);
FidEx0Vec = zeros(numel(NhV),1); EneEx0Vec = zeros(numel(NhV),1); Ene0ErrExVec = zeros(numel(NhV),1);

FidInVec = zeros(numel(NhV),1); EneInVec = zeros(numel(NhV),1); EneErrInVec = zeros(numel(NhV),1);
FidIn0Vec = zeros(numel(NhV),1); EneIn0Vec = zeros(numel(NhV),1); Ene0ErrInVec = zeros(numel(NhV),1);

gs = gs_ms; hmt = H_ms; Sz0Inds = (sum(basis,2)==0);
gs0 = gs(Sz0Inds); gs0 = gs0 / sqrt(sum(abs(gs0).^2)); hmt0 = hmt(Sz0Inds,Sz0Inds);

Param0 = AnsatzObj.Modifier{1}.ParamList; W0 = AnsatzObj.Modifier{1}.W; 
a0 = AnsatzObj.Modifier{1}.a; b0 = AnsatzObj.Modifier{1}.b;

HilbertObj = Spin(N,1/2,0); GraphObj = HypCub(N,1,1,1);
ModParams.nmag = 0; ModParams.nphs = 0; ModParams.a = 0; ModParams.b = 0; ModParams.W = 0;
ModParams.Nh = Nh-1; NQSEx = NQS(HilbertObj,GraphObj,ModParams,1);
ModParams.Nh = 1; NQSIn = NQS(HilbertObj,GraphObj,ModParams,1);

PsiFull = AnsatzObj.Modifier{1}.PsiGenerate(basis);
PsiFull0 = PsiFull(Sz0Inds); PsiFull0 = PsiFull0 / sqrt(sum(abs(PsiFull0).^2));
FidFull = abs(PsiFull'*gs)^2; FidFull0  = abs(PsiFull0'*gs0)^2;

figure(1); plot(NhV,FidFull*ones(numel(NhV),1),'--'); hold on;
figure(2); plot(NhV,FidFull0*ones(numel(NhV),1),'--'); hold on;
for n = NhV
    ParamEx = [a0; b0(NhV~=n); reshape(W0(NhV~=n,:),N*(Nh-1),1)];
    ParamIn = [a0; b0(n); W0(n,:).'];
    NQSEx = NQSEx.ParamLoad(ParamEx); PsiEx = NQSEx.PsiGenerate(basis);
    NQSIn = NQSIn.ParamLoad(ParamIn); PsiIn = NQSIn.PsiGenerate(basis);
    PsiEx0 = PsiEx(Sz0Inds); PsiEx0 = PsiEx0 / sqrt(sum(abs(PsiEx0).^2));
    PsiIn0 = PsiIn(Sz0Inds); PsiIn0 = PsiIn0 / sqrt(sum(abs(PsiIn0).^2));
    FidExVec(n) = abs(PsiEx'*gs)^2; FidEx0Vec(n) = abs(PsiEx0'*gs0)^2;
    FidInVec(n) = abs(PsiIn'*gs)^2; FidIn0Vec(n) = abs(PsiIn0'*gs0)^2;
    EnEx = PsiEx'*hmt*PsiEx; EnEx0 = PsiEx0'*hmt0*PsiEx0;
    EnIn = PsiIn'*hmt*PsiIn; EnIn0 = PsiIn0'*hmt0*PsiIn0;
    EneErrEx = abs((EnEx - en)/en); EneErrEx0 = abs((EnEx0 - en)/en);
    EneErrIn = abs((EnIn - en)/en); EneErrIn0 = abs((EnIn0 - en)/en);
    EneExVec(n) = EnEx; EneEx0Vec(n) = EnEx0; 
    EneInVec(n) = EnIn; EneIn0Vec(n) = EnIn0;
    EneErrExVec(n) = EneErrEx; Ene0ErrExVec(n) = EneErrEx0;
    EneErrInVec(n) = EneErrIn; Ene0ErrInVec(n) = EneErrIn0;
end
figure(1); plot(NhV,FidExVec,'-o'); plot(NhV,FidInVec,'-x'); 
legend('Full NQS','Excluding hidden unit','Only 1 hidden unit');
figure(2); plot(NhV,FidEx0Vec,'-o'); plot(NhV,FidIn0Vec,'-x');
legend('Full NQS','Excluding hidden unit','Only 1 hidden unit');