% Script for constructing and optimising Ansatz objects with Jastrow
% Modifiers and Bose condensate References for the Bose Hubbard model on a
% 2D lattice with periodic boundary conditions.

%% Specify the basic information for the lattice.
L = 4; N = L^2; Dim = [L L];

%% Specify Graph object.
LVecs = [1 0; 0 1]; SFlag = 1;
% SFlag determines whether to generate all possible separations in GraphObj.BondMap.
GraphObj = HypCub(Dim,'PBC',LVecs,SFlag);
% GraphObj describes a 2D L x L lattice with N sites and periodic boundary
% conditions, and contains a lookup table for nearest neighbours in x and
% y, as well as for multiple nearest neighbour separations.

NVec = [N-1 N N+1]; NStr = {' N-1 ',' N ',' N+1 '};

for nb = 2 % 1:numel(NVec)
    
    %% Specify Hilbert object.
    Nb = NVec(nb); Nmax = 4;
    HilbertObj = Bose(N,Nb,Nmax);
    
    %% Specify initial Hamiltonian object.
    dRT = [1 0; 0 1]; dRU = [0 0];
    % Set up the Graphs the Hamiltonian object will assign to its Operators.
    HGraphs = [HypCub(Dim,'PBC',dRT,1); HypCub(Dim,'PBC',dRU,0)];
    % Specify the Hamiltonian energy terms e.g. hopping, interaction, exchange.
    t = -1; U = [0 2 4 6 8 9 10 11 12 14 16 18 20 22 24 28 32 40 48 64 80 96 120 150];
    
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
    SamplerObj1 = Sampler(HilbertObj,HamiltonianObj,{});
    
    SamplerObj1 = SetNsamp(SamplerObj1,8000);
    
    SamplerObj2 = SetNsamp(SamplerObj1,16000);
    
    %% Initialise Stochastic Reconfig Optimiser.
    Npass1 = 300; Npass2 = 200; Ncore = 8; LR1 = 0.1; LR2 = 0.01;
    
    SR1 = StochasticReconfig(Npass1,Ncore);
    SR2 = StochasticReconfig(Npass2,Ncore);
    
    SR1 = SetExtraSamples(SR1,12000);
    SR2 = SetExtraSamples(SR2,24000);
    
    SR1 = SetSRTolerances(SR1,5+2*U(1),1e-6);
    SR2 = SetSRTolerances(SR2,5+2*U(1),1e-6);
    
    SR1 = SetRegularisation(SR1,1e4,1e-3,0.9);
    SR2 = SetRegularisation(SR2,1e-3,1e-3,1);
    
    SR1 = SetLearnRate(SR1,LR1);
    SR2 = SetLearnRate(SR2,LR2);
    
    AnsStr0 = 'BECR-Jast-NNMB SeqOpt';
    AnsStr = 'BECR-Jast-NNMB FT'; % Ansatz name string for easy identification later.
    
    urange{1} = [];
    
    urange{2} = 1:numel(U);
    
    urange{3} = [];
    
    % mkdir([AnsStr NStr{nb} num2str(L) 'x' num2str(L)]);
    
    for u = urange{nb}
        load(['BHM 2D U ' num2str(U(u)) NStr{nb} num2str(N) ' ' AnsStr0 ' Logs.mat'],'AnsatzObj');
        if u > 1
            % Alter Hamiltonian parameter corresponding to interaction.
            HamiltonianObj.HParams(2) = U(u)/2;
            SamplerObj1 = SetHamiltonian(SamplerObj1,HamiltonianObj);
            SamplerObj2 = SetHamiltonian(SamplerObj2,HamiltonianObj);
            % Reset energy tolerance in SR to avoid getting stuck early in
            % optimisation.
            SR1 = SetSRTolerances(SR1,5+2*U(u),1e-6);
            SR2 = SetSRTolerances(SR2,5+2*U(u),1e-6);
        end
        tic;
        % Perform optimisation of AnsatzObj with SR.
        EnIter = cell(2,1);
        [AnsatzObj,EnIter{1}] = SR1.Optimise(SamplerObj1,AnsatzObj);
        [AnsatzObj,EnIter{2}] = SR2.Optimise(SamplerObj2,AnsatzObj);
        EnIter = cell2mat(EnIter);
        RunTime = toc;
        save(['BHM 2D U ' num2str(U(u)) NStr{nb} num2str(N) ' ' AnsStr ' Logs.mat'],...
            'AnsatzObj','EnIter','RunTime');
        % Save AnsatzObj and run details for later analysis.
    end
end
