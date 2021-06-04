% --- Continuous time Green's Function Monte Carlo sampling function ---

function [EvalAvgP,WAvgP] = GFMCSampleCT(GFMCObj,AnsatzObj)

% Enacts importance sampling CTGFMC using provided Ansatz object for the
% guiding wavefunction.
Ncore = GFMCObj.Ncore; Nwalk = GFMCObj.Nwalk; Nsamp = GFMCObj.Nsamp;
Nequil = GFMCObj.Nequil; Pmax = GFMCObj.Pmax; Tbranch = GFMCObj.Tbranch; Teq = GFMCObj.Tequil;
HamiltonianObj = GFMCObj.Hamiltonian; HilbertObj = HamiltonianObj.Operator{1}.Hilbert;

Operators = GFMCObj.Operators; OpNum = numel(Operators);

EvalAvgP = cell(Nsamp,OpNum+1); WAvgP = zeros(Nsamp+Pmax,1);

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
                GFMCStepCT(HamiltonianObj,AnsW{w},CfgW{w},Teq,Weights(w),{});
        end
    else
        for w = 1:Nwalk
            [~,~,AnsW{w},CfgW{w},Weights(w)] = ...
                GFMCStepCT(HamiltonianObj,AnsW{w},CfgW{w},Teq,Weights(w),{});
        end
    end
    [CfgW,AnsW] = GFMCReconfig(CfgW,AnsW,Weights);
end

for n = 1:(Nsamp+Pmax)
    EneLocVec = zeros(Nwalk,1); OpValCell = cell(Nwalk,OpNum); Weights = ones(Nwalk,1);
    if Ncore > 1
        parfor (w=1:Nwalk,Ncore)
            [EneLocVec(w),OpValCell(w,:),AnsW{w},CfgW{w},Weights(w)] = ...
                GFMCStepCT(HamiltonianObj,AnsW{w},CfgW{w},Tbranch,Weights(w),Operators);
        end
    else
        for w = 1:Nwalk
            [EneLocVec(w),OpValCell(w,:),AnsW{w},CfgW{w},Weights(w)] = ...
                GFMCStepCT(HamiltonianObj,AnsW{w},CfgW{w},Tbranch,Weights(w),Operators);
        end
    end
    [EnAvg,OpAvg,CfgW,AnsW,WAvg,~] = GFMCReconfigEval(EneLocVec,OpValCell,CfgW,AnsW,Weights);
    WAvgP(n) = WAvg;
    if n > Pmax
        EvalAvgP{n-Pmax,1} = EnAvg;
        if OpNum > 0
            EvalAvgP(n-Pmax,2:end) = OpAvg;
        end
    end
end

if Ncore > 1
    delete(Pool);
end

end
