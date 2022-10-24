SaveStr = ['BHM 2D U ' num2str(U(u)) NStr num2str(N) ' ' AnsStr ' Logs.mat'];
% Alter Hamiltonian parameter corresponding to interaction.
HamiltonianObj.HParams(1) = t/U(u);
EvalHamiltonian.HParams(2) = U(u)/2;
SamplerObj1 = SetHamiltonian(SamplerObj1,HamiltonianObj);
SamplerObj2 = SetHamiltonian(SamplerObj2,HamiltonianObj);
EvalSampler = SetHamiltonian(EvalSampler,EvalHamiltonian);
% Reset energy tolerance in SR to avoid getting stuck early in
% optimisation.
SR1 = SetSRTolerances(SR1,5/U(u)+2,1e-6); SR1 = SetEnergyTolerances(SR1,dEV1/U(u),0.05);
SR1 = SetRegularisation(SR1,RegMax1(u),RegMin1(u),RegFac1(u)); SR1 = SetLearnRate(SR1,LearnRate1(u));
SR2 = SetSRTolerances(SR2,5/U(u)+2,1e-6); SR2 = SetEnergyTolerances(SR2,dEV2/U(u),0.05);
SR2 = SetRegularisation(SR2,RegMax2(u),RegMin2(u),RegFac2(u)); SR2 = SetLearnRate(SR2,LearnRate2(u));
tic;
% Perform optimisation of AnsatzObj with SR.
EnIter = cell(2,1);
[AnsatzObj,EnIter{1}] = SR1.Optimise(SamplerObj1,AnsatzObj);
[AnsatzObj,EnIter{2}] = SR2.Optimise(SamplerObj2,AnsatzObj);
EnIter = cell2mat(EnIter); Params = AnsatzObj.ParamList; RunTime = toc;
save([DirStr '/' SaveStr],'AnsatzObj','EnIter','RunTime','Params');
% Perform a sampling run to evaluate observables.
tic;
[EnAvg,EnSamp,EvalAvg,~] = MultiChainSample(EvalSampler,AnsatzObj,Ncore);
EvalTime = toc; RunDate = date;
EneGS = EnAvg/N; VarN = EvalAvg{1}; NiNj = reshape(EvalAvg{2},L,L); 
EnSamp = EnSamp(:)/N;
for b = 1:Nbin
    EneBin(b) = mean(EnSamp((1:Lbin)+(b-1)*Lbin));
end
VarE = mean((EneBin - mean(EneBin)).^2)/Nbin;
disp(['EnAvg: ' num2str(EnAvg) '  |  Bin average: ' num2str(mean(EneBin)) ]);
DbHl = reshape(EvalAvg{3},L,L); OcFr = reshape(EvalAvg{4},1,Nmax+1); BiBj = EvalAvg{5}; 
DiDj = reshape(EvalAvg{6},L,L); HiHj = reshape(EvalAvg{7},L,L); 
save([DirStr '/' SaveStr],'AnsatzObj','RunTime','EnIter','EneGS','VarN','NiNj','VarE','DbHl','OcFr','BiBj','DiDj','HiHj','Params','EvalTime','RunDate');
% Save AnsatzObj and run details for later analysis.
