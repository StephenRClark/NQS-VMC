% Script for constructing and optimising Ansatz objects with number-hidden
% NQS Modifiers and Bose condensate References for the Bose Hubbard model
% on a 1D lattice with periodic boundary conditions.

L = 6; N = L^2; Dim = [L L]; Nmax = 4; HilbertObj = Bose(N,N,Nmax);
TGraph = HypCub(Dim,[1 1],eye(2),1); HopOp = Operator2S(HilbertObj,TGraph,@BpBm_OpMatEls); 
UGraph = HypCub(Dim,[1 1],[0 0],1); IntOp = OperatorDg(HilbertObj,UGraph,@NiNj_Bose_CfgVal);

EnSqOp = OperatorDg(HilbertObj,UGraph,@TE_EnLocSqCorr);
VarNOp = OperatorDg(HilbertObj,UGraph,@VarN_Bose_CfgVal);
NiNjOp = OperatorDg(HilbertObj,TGraph,@NiNj_Bose_CfgVal);
DbHlOp = OperatorDg(HilbertObj,TGraph,@DbHl_Bose_CfgVal);
OcFrOp = OperatorDg(HilbertObj,UGraph,@OccFrac_Bose_CfgVal);
BiBjOp = Operator2S(HilbertObj,TGraph,@BpBm_OpMatEls);

U = [1 4 8 12 14 15 16 17 18 19 20 21 22 24 28 32 40 48]; t = -1;

Rmax1 = 1e5*ones(numel(U),1);
Rmin1 = 1e1*ones(numel(U),1);
Rfac1 = 0.953;
LR1 = 0.02*ones(numel(U),1);
Rmax2 = 1e3*ones(numel(U),1);
Rmin2 = 1e-1*ones(numel(U),1);
Rfac2 = 0.982;
LR2 = 0.05*ones(numel(U),1);
Rmax3 = 1e1*ones(numel(U),1);
Rmin3 = 1e-3*ones(numel(U),1);
Rfac3 = 0.953;
LR3 = 0.01*ones(numel(U),1);
ETol = 16*N; EvalSamples = 40000;
Nsamp1 = 6400; Nsamp2 = 8000; Nsamp3 = 9600;
ExtraSamp1 = 1600; ExtraSamp2 = 2000; ExtraSamp3 = 2400;
Npass1 = 200; Npass2 = 600; Npass3 = 200;
% Nsamp1 = 12800; Nsamp2 = 16000; Nsamp3 = 19200;
% ExtraSamp1 = 3200; ExtraSamp2 = 4000; ExtraSamp3 = 4800;
% Npass1 = 100; Npass2 = 300; Npass3 = 100;
   
dEV1 = 0.2; dEV2 = 0.1; dEV3 = 0.05;

urange = 1:numel(U);

for u = urange
    NumSampler = 3;
    HamiltonianObj = Hamiltonian({HopOp;IntOp},[t/U(u), 1/2]);
    EvalHamiltonian = Hamiltonian({HopOp;IntOp},[t, U(u)/2]);
    EvalSampler = Sampler(HilbertObj,EvalHamiltonian,{});
    EvalSampler = EvalSampler.SetNsamp(EvalSamples); EvalSampler = EvalSampler.SetNblock(N);
    Samp1 = Sampler(HilbertObj,HamiltonianObj,{}); 
    Samp1 = Samp1.SetNsamp(Nsamp1); Samp1 = Samp1.SetNblock(N);
    Samp2 = Sampler(HilbertObj,HamiltonianObj,{}); 
    Samp2 = Samp2.SetNsamp(Nsamp2); Samp2 = Samp2.SetNblock(N);
    Samp3 = Sampler(HilbertObj,HamiltonianObj,{}); 
    Samp3 = Samp3.SetNsamp(Nsamp3); Samp3 = Samp3.SetNblock(N);
    SR1 = StochasticReconfig(Npass1,1); SR1 = SR1.SetExtraSamples(ExtraSamp1);
    SR1 = SR1.SetRegularisation(Rmax1(u),Rmin1(u),Rfac1);
    SR1 = SR1.SetSRTolerances(ETol,1e-6); SR1 = SR1.SetLearnRate(LR1(u));
    SR2 = StochasticReconfig(Npass2,1); SR2 = SR2.SetExtraSamples(ExtraSamp2);
    SR2 = SR2.SetRegularisation(Rmax2(u),Rmin2(u),Rfac2);
    SR2 = SR2.SetSRTolerances(ETol,1e-6); SR2 = SR2.SetLearnRate(LR2(u));
    SR3 = StochasticReconfig(Npass3,1); SR3 = SR3.SetExtraSamples(ExtraSamp3);
    SR3 = SR3.SetRegularisation(Rmax3(u),Rmin3(u),Rfac3);
    SR3 = SR3.SetSRTolerances(ETol,1e-6); SR3 = SR3.SetLearnRate(LR3(u));
    Samp1Prop = Samp1.PropertyList; Samp2Prop = Samp2.PropertyList; Samp3Prop = Samp3.PropertyList;
    SR1Prop = SR1.PropertyList; SR2Prop = SR2.PropertyList; SR3Prop = SR3.PropertyList; 
    HilbertProp = HilbertObj.PropertyList; EvalSampProp = EvalSampler.PropertyList;
    EnSqProp = EnSqOp.PropertyList; VarNProp = VarNOp.PropertyList; NiNjProp = NiNjOp.PropertyList;
    DbHlProp = DbHlOp.PropertyList; OcFrProp = OcFrOp.PropertyList; BiBjProp = BiBjOp.PropertyList;
    save(['BHM 2D U ' num2str(U(u)) ' N ' num2str(N) ' optimisation set 4 setup.mat'],...
        'HilbertProp','NumSampler','SR1Prop','SR2Prop','SR3Prop','Samp1Prop','Samp2Prop','Samp3Prop',...
        'EvalSampProp','EnSqProp','VarNProp','NiNjProp','DbHlProp','OcFrProp','BiBjProp');
end

save(['BHM 2D N ' num2str(N) ' optimisation set 4 parameters.mat'],'U','N',...
    'HilbertObj','Rmax1','Rmin1','Rfac1','Rmax2','Rmin2','Rfac2','Rmax3','Rmin3','Rfac3',...
    'ETol','Nsamp1','Nsamp2','Nsamp3','Npass1','Npass2','Npass3','LR1','LR2','LR3',...
    'EvalSamples','ExtraSamp1','ExtraSamp2','ExtraSamp3','dEV1','dEV2','dEV3');