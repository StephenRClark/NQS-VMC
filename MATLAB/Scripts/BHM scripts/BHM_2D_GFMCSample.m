% Script for sampling pre-optimised Ansatz objects with the Bose Hubbard
% model and any specified Operators.

%% Specify the basic information for the lattice.
L = 3; Dim = [L L]; N = prod(Dim);

%% Specify Graph object.
LVecs = eye(2); SFlag = 1; Bound = [1 1];
% SFlag determines whether to generate all possible separations in GraphObj.BondMap.
GraphObj = HypCub(Dim,Bound,LVecs,SFlag);
% GraphObj describes a 2D lattice with N sites and periodic boundary
% conditions, and contains a lookup table for nearest neighbours, as well
% as for multiple nearest neighbour separations.

%% Specify Hilbert object.
Nb = N; Nmax = 3;
HilbertObj = Bose(N,Nb,Nmax);

%% Specify initial Hamiltonian object.
dRT = eye(2); dRU = [0 0];
% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,Bound,dRT,1); HypCub(Dim,Bound,dRU,0)];
% Specify the Hamiltonian energy terms e.g. hopping, interaction, exchange.
t = -1; U = [1 4 8 12 14 15 16 17 18 19 20 21 22 24 28 32 40 48];
% Specify the Operators that will be used in the Hamiltonian.
HopOp = Operator2S(HilbertObj,HGraphs(1),@BpBm_OpMatEls);
% HopOp describes the bosonic hopping term, and uses the nearest neighbour
% lookup table provided by HGraphs(1).
IntOp = OperatorDg(HilbertObj,HGraphs(2),@NiNj_Bose_CfgVal);
% IntOp describes the density-density interaction term, and uses the
% entirely local lookup table provided by HGraphs(2)
Operators = {HopOp; IntOp}; HParams = [t; U(1)/2];

HamiltonianObj = Hamiltonian(Operators,HParams);

%% Initialise Sampler.

AnsStr = 'BECR-Jast'; u = 4;

load(['BHM 2D U ' num2str(U(u)) ' N ' num2str(N) ' ' AnsStr ' Logs.mat'],'AnsatzObj');

% HamiltonianObj.HParams(1) = t/U(u);
HamiltonianObj.HParams(2) = U(u)/2; 

Ncore = 1; Nwalk = 50; Pmax = Nwalk; Tbranch = 2; Nsamp = 250/Tbranch;  Nequil = 1000/Tbranch;

GFMCObj = GFMCSampler(Ncore,HamiltonianObj,{});

GFMCObj = GFMCObj.SetNwalk(Nwalk,Ncore); GFMCObj = GFMCObj.SetNequil(Nequil); GFMCObj.SetTequil(Tbranch);
GFMCObj = GFMCObj.SetNsamp(Nsamp); GFMCObj = GFMCObj.SetProjection(Pmax); GFMCObj = GFMCObj.SetTbranch(Tbranch);

[EvalAvg,WAvg] = GFMCObj.GFMCSampleCT(AnsatzObj);

WProj = ones(Nsamp,Pmax); EnAvg = cell2mat(EvalAvg);

load(['BHM 2D N ' num2str(N) ' U ' num2str(U(u)) ' exact ground state.mat'],'gs_en');

for p = 1:Pmax
    for q = p:Pmax
        WProj(:,q) = WProj(:,q) .* WAvg((1:Nsamp)+Pmax+1-p);
    end
end
EnSum = WProj .* EnAvg; EnProj = sum(EnSum,1)./sum(WProj,1);
% plot(0:Pmax,gs_en*ones(Pmax+1,1),'--'); hold on; 
plot(1:Pmax,EnProj,'-o'); 