% --- Continuous time Green's Function Monte Carlo Markov chain generating function ---

function [EnAvgP,WAvgP,CfgsP,RcfgIndsP] = GFMCChain(GFMCObj,AnsatzObj)

% Enacts importance sampling CTGFMC using provided Ansatz object for the
% guiding wavefunction.
Ncore = GFMCObj.Ncore; Nwalk = GFMCObj.Nwalk; Nsamp = GFMCObj.Nsamp;
Nequil = GFMCObj.Nequil; Pmax = GFMCObj.Pmax; T = GFMCObj.T;
HilbertObj = AnsatzObj.Hilbert; HamiltonianObj = GFMCObj.Hamiltonian;

EnAvgP = zeros(Nsamp+Pmax,1); WAvgP = zeros(Nsamp+Pmax,1);
RcfgIndsP = zeros(Nsamp+Pmax,Nwalk); CfgsP = cell(Nsamp+Pmax,Nwalk);

% Set up walker-specific properties.
CfgW = cell(Nwalk,1); AnsW = cell(Nwalk,1);
for w = 1:Nwalk % Generate Cfgs and individual copies of Ansatz per walker.
    CfgW{w} = HilbertObj.RandomCfg; AnsW{w} = PrepPsi(AnsatzObj,CfgW{w});
end

if Ncore > 1
    Pool = parpool(Ncore); % Save for deletion at the end.
end

for n = 1:Nequil
    Weights = ones(Nwalk,1);
    if Ncore > 1
        parfor (w=1:Nwalk,Ncore)
            [~,~,AnsW{w},CfgW{w},Weights(w)] = ...
                GFMCStepCT(HamiltonianObj,AnsW{w},CfgW{w},T,Weights(w),{});
        end
    else
        for w = 1:Nwalk
            [~,~,AnsW{w},CfgW{w},Weights(w)] = ...
                GFMCStepCT(HamiltonianObj,AnsW{w},CfgW{w},T,Weights(w),{});
            
        end
    end
    % Enact reconfiguration amongst walkers.
    [CfgW,AnsW] = GFMCReconfig(CfgW,AnsW,Weights);
end

for n = 1:(Nsamp + Pmax)
    EneLocVec = zeros(Nwalk,1); Weights = ones(Nwalk,1);
    if Ncore > 1
        parfor (w=1:Nwalk,Ncore)
            [EneLocVec(w),~,AnsW{w},CfgW{w},Weights(w)] = ...
                GFMCStepCT(HamiltonianObj,AnsW{w},CfgW{w},T,Weights(w),{});
        end
    else
        for w = 1:Nwalk
            [EneLocVec(w),~,AnsW{w},CfgW{w},Weights(w)] = ...
                GFMCStepCT(HamiltonianObj,AnsW{w},CfgW{w},T,Weights(w),{});
        end
    end
    % Enact reconfiguration amongst walkers and include indices.
    [EnAvg,~,CfgW,AnsW,WAvg,Inds] = GFMCReconfigEval(EneLocVec,{},CfgW,AnsW,Weights);
    EnAvgP(n) = EnAvg; WAvgP(n) = WAvg;
    RcfgIndsP(n,:) = Inds; CfgsP(n,:) = CfgW;
end

if Ncore > 1
    delete(Pool);
end

end
