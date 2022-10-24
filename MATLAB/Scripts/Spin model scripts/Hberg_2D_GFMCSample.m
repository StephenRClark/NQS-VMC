% Script for sampling pre-optimised Ansatz objects with the Heisenberg
% model and any specified Operators.

%% Specify the basic information for the lattice.
L = 4; N = L^2; Dim = [L L];

%% Specify Graph object.
LVecs = eye(numel(Dim)); SFlag = 1;
% SFlag determines whether to generate all possible separations in GraphObj.BondMap.
GraphObj = HypCub(Dim,'PBC',LVecs,SFlag);
% GraphObj describes a 2D lattice with N sites and periodic boundary
% conditions, and contains a lookup table for nearest neighbours, as well
% as for multiple nearest neighbour separations.
for i = 1:L
    for j = 1:L
        GraphObj.SLInds(i+(j-1)*L) = 1 + mod(i+j,2);
    end
end

%% Specify Hilbert object.
S = 1/2; Sector = 0;
HilbertObj = Spin(N,S,Sector);

%% Specify initial Hamiltonian object.

J = 1;

SiSjOp = Operator2S(HilbertObj,GraphObj,@SiSj_GT_OpMatEls);

Operators = {SiSjOp}; HParams = J;

HamiltonianObj = Hamiltonian(Operators,HParams);

%% Initialise Sampler.

Lambda = 1000; % Projection 'reference energy'.

Nequil = 5000; Ncore = 4; Nsamp = 1000; % Basic sampler parameters.

Nproj = 10; Nwalk = 20; Nbranch = 5; Tbranch = 0.1; % GFMC specific parameters.

GFMCObj = GFMCSampler(Ncore,HamiltonianObj,{});

GFMCObj = SetNwalk(GFMCObj,Nwalk,Ncore); GFMCObj = SetProjection(GFMCObj,Nproj);
GFMCObj = SetNsamp(GFMCObj,Nsamp); GFMCObj = SetBranch(GFMCObj,Nbranch,Tbranch);
GFMCObj = SetNequil(GFMCObj,Nequil); GFMCObj = SetRefEnergy(GFMCObj,Lambda);

%% Select Ansatz for guiding wavefunction.

AnsStr = 'Plus-NQSTI Alpha 1';

load(['Heisenberg 2D N ' num2str(N) ' ' AnsStr ' Logs.mat'],'AnsatzObj');
% Load exact results for comparison.
% load(['Heisenberg 2D N ' num2str(N) ' exact ground state.mat']);%,'gs_en');

% Begin GFMC sampling.
[EnAvgP,WAvgP] = GFMCObj.GFMCSample(AnsatzObj);
% Projection index analysis setup.
EnAvgP = cell2mat(EnAvgP);
ProjWeights = zeros(Nsamp,Nproj);
for p = 1:Nproj
    for n = 1:Nsamp
        ProjWeights(n,p) = prod(WAvgP(n + ((Nproj-p+1):Nproj)));
    end
end
E0P = EnAvgP((1:Nsamp)+Nproj)'*ProjWeights ./ sum(ProjWeights,1);
% Plot calculated energies using different projection indices.
figure(1); plot(E0P/N,'-o'); % hold on; plot(gs_en*ones(Nproj,1),'-x');
