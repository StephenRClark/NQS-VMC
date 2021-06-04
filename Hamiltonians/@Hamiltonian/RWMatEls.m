% --- Weighted operator matrix element sampling function ---

function [CfgP, RWOpCells, CfgIndsP] = RWMatEls(HamiltonianObj,AnsatzObj,Cfg,InitCells,CfgInds)
% This function takes some input configurations and associated values (such
% as operator local values) and outputs the set of configurations linked by
% the input Hamiltonian as well as operator local values weighted by
% PsiRatio. Each CfgP has an associated RWOpCells cell array containing
% accrued operator local values which are multiplied later.

% Hamiltonian objects contain Operators which themselves have RWMatEls
% methods, so just invoke those and multiply by the relevant HParams.

HInd = size(InitCells,2)+1; % Keep track of which values of RWOpCells are related to H.

HOperators = HamiltonianObj.Operator; HParams = HamiltonianObj.HParams;

HCells = cell(numel(HParams),3); CCount = zeros(numel(HParams),1);

for h = 1:numel(HParams)
    if strcmp(class(HOperators{h}),'Operator2S')
        Ind = zeros(2,1); % Ensure correct formatting of Ind for RWMatEls.
    else
        Ind = 0;
    end
    [HCfg, HRWOpCells, HCfgIndsP] = RWMatEls(HOperators{h},AnsatzObj,Cfg,InitCells,CfgInds,Ind);
    CCount(h) = numel(HCfg);
    for c = 1:CCount(h) % Multiply by relevant Hamiltonian parameter.
        HRWOpCells{c,HInd} = HRWOpCells{c,HInd}.*HParams(h);
    end
    HCells(h,:) = {HCfg, HRWOpCells, HCfgIndsP};
end

% Initialise structure array for CfgP and arrays for CfgIndsP and
% RWOpCells.
CfgP(sum(CCount)) = Cfg; CfgP = reshape(CfgP,numel(CfgP),1);
RWOpCells = cell(sum(CCount),HInd); CfgIndsP = zeros(sum(CCount),size(CfgInds,2));

% Recombine contributions from each Operator.
for h = 1:numel(HParams)
    CInds = (1:CCount(h)) + sum(CCount(1:(h-1)));
    CfgP(CInds) = HCells{h,1};
    RWOpCells(CInds,:) = HCells{h,2};
    CfgIndsP(CInds,:) = HCells{h,3};
end
end