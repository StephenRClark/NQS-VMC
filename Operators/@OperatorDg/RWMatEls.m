% --- Single configuration correlation profile sampling function ---

function [CfgP, RWOpCells, CfgInds] = RWMatEls(OperatorObj,~,Cfg,InitCells,CfgInds,Ind)
% This function takes some input configurations and associated values (such
% as operator local values) and outputs the set of configurations linked by
% the input Operator as well as operator local values weighted by PsiRatio.
% Each CfgP has an associated RWOpCells cell array containing accrued
% operator local values which are multiplied later.

% For a diagonal operator, no configurations change, and the input values
% are simply multiplied by the configuration value.

HilbertObj = OperatorObj.Hilbert; GraphObj = OperatorObj.Graph;
CfgP = Cfg; TempCells = cell(numel(CfgP),1);
for c = 1:numel(CfgP)
    % Some operators output a vector naturally - need to distinguish
    % these cases from the size of CfgVal.
    CfgValT = OperatorObj.CfgVal(HilbertObj,CfgP(c),0,0,GraphObj.Bonds);
    if Ind ~= 0 && numel(CfgValT) == 1 
        % Scalar outputs associated with a set distance two-site operator -
        % need to supply with all available separations to obtain full set
        % of relevant values.
        CfgVal = zeros(numel(GraphObj.BondMap),1);
        for b = 1:numel(GraphObj.BondMap)
            CfgVal(b) = OperatorObj.CfgVal(HilbertObj,CfgP(c),0,0,GraphObj.BondMap{b});
        end
    else
        CfgVal = reshape(CfgValT,numel(CfgValT),1);
    end
    if Ind > 1
        CfgVal = reshape(CfgVal,[ones(1,Ind-1) numel(CfgVal)]);
    end
    % Would normally be weighted by PsiRatio, but no configuration change.
    TempCells{c} = CfgVal;
end
RWOpCells = [InitCells TempCells];
% CfgInds tracks any particular sum order for any new configurations
% generated. For a diagonal operator, no new configurations are generated,
% so the sum order is unchanged.
end