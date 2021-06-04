% --- Continuous time Green's Function Monte Carlo sampling function ---

function [EvalAvgP,WAvgP] = GFMCSampleDT(GFMCObj,AnsatzObj)

% Enacts importance sampling CTGFMC using provided Ansatz object for the
% guiding wavefunction.
Ncore = GFMCObj.Ncore; Nwalk = GFMCObj.Nwalk; Nsamp = GFMCObj.Nsamp;
Nequil = GFMCObj.Nequil; Pmax = GFMCObj.Pmax; Lambda = GFMCObj.Lambda; Nbranch = GFMCObj.Nbranch;
HilbertObj = AnsatzObj.Hilbert; HamiltonianObj = GFMCObj.Hamiltonian;

Operators = GFMCObj.Operators; OpNum = numel(Operators);

EvalAvgP = cell(Nsamp+Pmax,OpNum+1); WAvgP = zeros(Nsamp+Pmax,1);

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
                GFMCStep(HamiltonianObj,AnsW{w},CfgW{w},Lambda,Weights(w),{});
        end
    else
        for w = 1:Nwalk
            [~,~,AnsW{w},CfgW{w},Weights(w)] = ...
                GFMCStep(HamiltonianObj,AnsW{w},CfgW{w},Lambda,Weights(w),{});
        end
    end
    [CfgW,AnsW] = GFMCReconfig(CfgW,AnsW,Weights);
end

for n = 1:(Nsamp+Pmax)
    EneLocVec = zeros(Nwalk,1); OpValCell = cell(Nwalk,OpNum); Weights = ones(Nwalk,1);
    for m = 1:Nbranch
        if Ncore > 1
            parfor (w=1:Nwalk,Ncore)
                [EneLocVec(w),OpValCell(w,:),AnsW{w},CfgW{w},Weights(w)] = ...
                    GFMCStep(HamiltonianObj,AnsW{w},CfgW{w},Lambda,Weights(w),Operators);
            end
        else
            for w = 1:Nwalk
                [EneLocVec(w),OpValCell(w,:),AnsW{w},CfgW{w},Weights(w)] = ...
                    GFMCStep(HamiltonianObj,AnsW{w},CfgW{w},Lambda,Weights(w),Operators);
            end
        end
    end
    [EnAvg,OpAvg,CfgW,AnsW,WAvg,~] = GFMCReconfigEval(EneLocVec,OpValCell,CfgW,AnsW,Weights);
    EvalAvgP{n,1} = EnAvg; WAvgP(n) = WAvg;
    if OpNum > 0
        EvalAvgP(n,2:end) = OpAvg;
    end
end

if Ncore > 1
    delete(Pool);
end

end