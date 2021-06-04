% --- Composite operator correlation profile sampling function ---

function [CorrSamp] = GraphSample(CompOperator,Cfg,~,~,AnsatzObj)
% This function evaluates the expectation values of a composite Operator
% over the principal Bonds of the associated Graph.
N = Cfg.N;
% Unpack useful information from CompOperator.
OpOrder = CompOperator.OpOrder; SubOperator = CompOperator.SubOperator;
CorrIndices = cell2mat(OpOrder(:,3)); IndMax = max(CorrIndices);
CfgInds = zeros(1,IndMax); % For bookkeeping in later sum over generated configurations.
LOpInds = find(strcmp(OpOrder(:,1),'L')); ROpInds = find(strcmp(OpOrder(:,1),'R'));
LOpOrder = cell2mat(OpOrder(LOpInds,2)); ROpOrder = cell2mat(OpOrder(ROpInds,2));

LCfg = Cfg; LCfgInds = CfgInds; LOpCells = {}; % Sample left operators first.
for l = 1:numel(LOpInds)
    LOpNum = LOpInds(LOpOrder == l);
    [LCfg, LOpCells, LCfgInds] = SubOperator{LOpNum}.RWMatEls(AnsatzObj,...
        LCfg,LOpCells,LCfgInds,OpOrder{LOpNum,3});
end

RCfg = Cfg; RCfgInds = CfgInds; ROpCells = {}; % Sample right operators second.
for l = 1:numel(ROpInds)
    ROpNum = ROpInds(ROpOrder == l);
    [RCfg, ROpCells, RCfgInds] = SubOperator{ROpNum}.RWMatEls(AnsatzObj,...
        RCfg,ROpCells,RCfgInds,OpOrder{ROpNum,3});
end

% Perform sum using cells, then element-wise multiply left and right
% contributions to finalise.

% Find relevant cell dimensions first.
LCellDims = max([ones(1,IndMax); N*ceil(max(LCfgInds)/N)]);
% Any zero entries in CfgInds replaced with 1.
LIndAncilla = cumprod([1, LCellDims(1:(IndMax-1))]);
% Need to linearise indices to properly access and modify cell contents.
LCells = cell(LCellDims);
for lc = 1:numel(LCfg)
    LCfgInds(lc,:) = max([LCfgInds(lc,:); ones(1,IndMax)]);
    LinLCfgInd = 1+sum((LCfgInds(lc,:)-1).*LIndAncilla); % Linearised index.    
    LCfgVal = 1; % Reassemble values associated with this configuration.
    for li = 1:size(LOpCells,2)
        LCfgVal = LCfgVal .* LOpCells{lc,li};
    end
    if isempty(LCells{LinLCfgInd})
        LCells{LinLCfgInd} = 0;
    end
    % Add contribution to LCells.
    LCells{LinLCfgInd} = LCells{LinLCfgInd} + LCfgVal;
end
for l = 1:numel(LCells)
    if isempty(LCells{l})
        LCells{l} = 0 * LCfgVal;
    end
end

% Repeat above procedure for right operators.
% Find relevant cell dimensions first.
RCellDims = max([ones(1,IndMax); N*ceil(max(RCfgInds)/N)]);
% Any zero entries in CfgInds replaced with 1.
RIndAncilla = cumprod([1, RCellDims(1:(IndMax-1))]);
% Need to linearise indices to properly access and modify cell contents.
RCells = cell(RCellDims);
for rc = 1:numel(RCfg)
    RCfgInds(rc,:) = max([RCfgInds(rc,:); ones(1,IndMax)]);
    LinRCfgInd = 1+sum((RCfgInds(rc,:)-1).*RIndAncilla); % Linearised index.
    if isempty(RCells{LinRCfgInd})
        RCells{LinRCfgInd} = 0;
    end
    RCfgVal = 1; % Reassemble values associated with this configuration.
    for ri = 1:size(ROpCells,2)
        RCfgVal = RCfgVal .* ROpCells{rc,ri};
    end
    % Add contribution to RCells.
    RCells{LinRCfgInd} = RCells{LinRCfgInd} + RCfgVal;
end
for r = 1:numel(RCells)
    if isempty(RCells{r})
        RCells{r} = 0 * RCfgVal;
    end
end        
% Final recombination involves conversion of LCells and RCells - if all
% goes to plan, should be a clean element-wise multiplication.
CorrSamp = conj(cell2mat(LCells)) .* cell2mat(RCells);
end