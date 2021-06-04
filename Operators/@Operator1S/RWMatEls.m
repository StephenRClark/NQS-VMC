% --- Single configuration correlation profile sampling function ---

function [CfgP, RWOpCells, CfgInds] = RWMatEls(OperatorObj,AnsatzObj,Cfg,InitCells,CfgInds,Ind)
% This function takes some input configurations and associated values (such
% as operator local values) and outputs the set of configurations linked by
% the input Operator as well as operator local values weighted by PsiRatio.
% Each CfgP has an associated RWOpCells cell array containing accrued
% operator local values which are multiplied later.

% For a non-diagonal one site operator, if any Ind is specified, each
% configuration associated with a particular site needs to be marked with
% CfgInds for bookkeeping in later summation.

HilbertObj = OperatorObj.Hilbert; 
CfgCells = cell(numel(Cfg),1); IndCells = mat2cell(CfgInds,ones(numel(Cfg),1));
TempCells = cell(numel(Cfg),1);
DCount = zeros(numel(Cfg),1);
for c = 1:numel(Cfg)
    if Ind ~= 0 % Free index of a one-site operator corresponds to individual site index.
        [Diff,OpMatEls] = OperatorObj.CorrMatEls(Cfg(c));
        % Individual OpMatEls are scalar.
        Diff = reshape(Diff,numel(Diff),1); IndList = zeros(numel(Diff),1);
        for d = 1:numel(Diff)
            % Weight each matrix element by the associated PsiRatio.
            PsiRatio = PsiRatio(AnsatzObj,Diff(d));
            if isnan(PsiRatio) || isinf(PsiRatio)
                PsiRatio = 0;
            end
            OpMatEls(d,:) = OpMatEls(d,:) * PsiRatio;
            % Each configuration needs to be marked according to the site
            % the operator acts on.
            IndList(d) = Diff(d).pos;
        end
        % All generated configurations here are associated with a
        % particular site. Need to log this in IndCells.
        % IndCells{c} is a 1xM vector - need to change entries at Ind.
        % Alternatively, if that entry is already filled, only keep Diffs
        % that have the same value of that entry.
        if IndCells{c}(Ind(1)) ~= 0
            Diff = Diff(IndList==IndCells{c}(Ind(1)));
            OpMatEls = OpMatEls(IndList==IndCells{c}(Ind(1)),:);
        end
        if numel(Diff) ~= 0
            DCount(c) = numel(Diff);
            IndCells{c} = ones(numel(Diff),1) .* IndCells{c};
            IndCells{c}(:,Ind(1)) = IndCells{c}(:,Ind(1)) + IndList;
            % Convert Diffs to Cfgs and store in CfgCells.
            CfgCells{c} = HilbertObj.Diff2Cfg(Diff,Cfg(c));
            % Store weighted matrix elements in TempCells for later
            % unpacking.
            if (numel(Ind) == 2) && (numel(OpMatEls) ~= 0)% If a second Ind is included, assume Ind(2) is specified.
                OpMatEls = reshape(OpMatEls,[numel(Diff) ones(1,Ind(2)-2) numel(OpMatEls(1,:))]);
            end
            TempCells{c} = OpMatEls;
        end
    else % Only need to consider principal bonds in Graph.
        [Diff,OpMatEls] = OperatorObj.CorrMatEls(Cfg(c));
        Diff = reshape(Diff,numel(Diff),1);
        for d = 1:numel(Diff)
            % Weight each matrix element by the associated PsiRatio.
            OpMatEls(d,:) = OpMatEls(d,:) * PsiRatio(AnsatzObj,Diff(d));
        end
        % Convert Diffs to Cfgs and store in CfgCells.
        CfgCells{c} = HilbertObj.Diff2Cfg(Diff,Cfg(c));
    end
end
% CfgInds tracks any particular sum order for any new configurations
% generated. Can regenerate the list from IndCells.
CfgInds = cell2mat(IndCells);
% Unpack all the generated configurations from CfgCells.
CfgP(sum(DCount)) = CfgCells{1}(1);
for c = 1:numel(Cfg)
    CfgP((1:DCount(c))+sum(DCount(1:(c-1)))) = CfgCells{c};
end
% Create storage for collection of accrued operator values.
RWOpCells = cell(numel(CfgP),size(InitCells,2)+1); rcount = 1;
for c = 1:numel(Cfg) % Loop over supplied configurations.
    for d = 1:size(TempCells{c},1) % Loop over generated configurations.
        CfgVal = TempCells{c}(d,:);
        if numel(Ind) == 2 % If a second Ind is included, assume Ind(2) is specified.
            CfgVal = reshape(CfgVal,[ones(1,Ind(2)-1) numel(CfgVal)]);
        end
        if isempty(InitCells)
            RWOpCells(rcount,:) = {CfgVal};
        else
            RWOpCells(rcount,:) = {InitCells{c,:}, CfgVal}; % Append new matrix elements to cell array.
        end
        rcount = rcount + 1;
    end
end

end