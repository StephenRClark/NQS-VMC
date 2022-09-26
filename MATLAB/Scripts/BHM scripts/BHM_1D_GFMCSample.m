% Script for sampling pre-optimised Ansatz objects with the Bose Hubbard
% model and any specified Operators.

%% Specify the basic information for the lattice.
N = 10; Dim = N;

%% Specify Graph object.
LVecs = 1; SFlag = 1;
% SFlag determines whether to generate all possible separations in GraphObj.BondMap.
GraphObj = HypCub(Dim,'PBC',LVecs,SFlag);
% GraphObj describes a 1D lattice with N sites and periodic boundary
% conditions, and contains a lookup table for nearest neighbours, as well
% as for multiple nearest neighbour separations.

%% Specify Hilbert object.
Nb = N; Nmax = 4;
HilbertObj = Bose(N,Nb,Nmax);

%% Specify initial Hamiltonian object.
dRT = 1; dRU =  0;
% Set up the Graphs the Hamiltonian object will assign to its Operators.
HGraphs = [HypCub(Dim,'PBC',dRT,1); HypCub(Dim,'PBC',dRU,0)];
% Specify the Hamiltonian energy terms e.g. hopping, interaction, exchange.
t = -1; U = 10;
% Reentrance U vector: [2 5/2 20/7 10/3 4 13/3 37/8 5 21/4 50/9 40/7 29/5 53/9 6 25/4 20/3 15/2 35/4 10 12 16 20 30 50 100]
% Exact test U vector: [0 1 2 3 4 5 5.5 6 6.5 7 8 9 10 12 16 24]

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

Nequil = 5000; Ncore = 4; Nsamp = 1000; % Basic sampler parameters.

Nproj = 100; Nwalk = 20; Tbranch = 0.1; % GFMC specific parameters.

GFMCObj = GFMCSampler(Ncore,HamiltonianObj,{});

GFMCObj = SetNwalk(GFMCObj,Nwalk,Ncore); GFMCObj = SetProjection(GFMCObj,Nproj);
GFMCObj = SetNsamp(GFMCObj,Nsamp); GFMCObj = SetTBranch(GFMCObj,Tbranch);
GFMCObj = SetNequil(GFMCObj,Nequil); 

%% Select Ansatz for guiding wavefunction.

AnsStr = 'BECR-Jast FT2'; NStr = ' N '; urange = 1:numel(U);

% mkdir([AnsStr NStr num2str(N) ' GFMC results']); EneGFMCU = zeros(numel(U),1);

for u = urange

HamiltonianObj.HParams(2) = U(u)/2; GFMCObj = SetHamiltonian(GFMCObj,HamiltonianObj);

load(['BHM 1D U ' num2str(U(u)) NStr num2str(N) ' ' AnsStr ' Logs.mat'],'AnsatzObj');

JastObj = AnsatzObj.Modifier{1}; BECRObj = AnsatzObj.Reference;

AnsatzObj = Ansatz(BECRObj,{JastObj},HilbertObj);

[EnAvg,WAvg] = GFMCObj.GFMCSample(AnsatzObj);
% Projection index analysis setup.
EnAvg = cell2mat(EnAvg);
ProjWeights = ones(Nsamp,Nproj+1);
ProjWeights(:,1) = WAvg((1:Nsamp)+Nproj);
for p = 1:Nproj
    for m = 1:p
        ProjWeights(:,p+1) = ProjWeights(:,p+1) .* WAvg((1:Nsamp)+Nproj-m);
    end
end
E0P = EnAvg((1:Nsamp)+Nproj)'*ProjWeights ./ sum(ProjWeights,1);
EneGFMCU(u) = E0P(1);
save(['BHM 1D U ' num2str(U(u)) NStr num2str(N) ' ' AnsStr ' GFMC Logs P100.mat'],...
    'EnAvg','WAvg','ProjWeights','E0P');
disp(['Sampling for U = ' num2str(U(u)) ' complete.']);
end

figure(1); plot(U(urange),EneGFMCU,'-o');
