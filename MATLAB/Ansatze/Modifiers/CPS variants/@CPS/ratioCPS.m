function [Ratio,Update] = ratioCPS(CPSTIObj,Diff)
    Diff.oldVal = zeros(Diff.num,1);
    Diff.newVal = zeros(Diff.num,1);
    for i = 1:Diff.num
        Diff.oldVal(i) = CPSTIObj.tempCfgVec(Diff.pos(i));
        Diff.newVal(i) = Diff.oldVal(i) + Diff.val(i);
    end
    thetaShift = CPSTIObj.Theta;
    for j = 1:CPSTIObj.Nh
        for k = 1:Diff.num
            for l = 1:2
                if Diff.oldVal(k) ~= 0
                    thetaShift((j - 1) * 2 + l) = thetaShift((j - 1) * 2 + l) - CPSTIObj.W((j - 1) * CPSObj.Nv * 4 + (Diff.pos(k) - 1) * 4 + (l - 1) * 2 + Diff.oldVal(k));
                end
                if Diff.newVal(k) ~= 0
                    thetaShift((j - 1) * 2 + l) = thetaShift((j - 1) * 2 + l) + CPSTIObj.W((j - 1) * CPSObj.Nv * 4 + (Diff.pos(k) - 1) * 4 + (l - 1) * 2 + Diff.newVal(k));
                end
            end
        end
    end
    Ratio = 0;
    for k = 1:Diff.num
        if Diff.oldVal(k) ~= 0
            Ratio = Ratio - a((Diff.pos(k) - 1) * 2 + Diff.oldVal(k));
        end
        if Diff.newVal(k) ~= 0
            Ratio = Ratio + a((Diff.pos(k) - 1) * 2 + Diff.newVal(k));
        end
    end
    Ratio = exp(Ratio);
    for j = 1:CPSTIObj.Nh
        Ratio = Ratio * ((1 + thetaShift((j - 1) * 2 + 1) + thetaShift(j * 2)) / (1 + CPSTIObj.Theta((j - 1) * 2 + 1) + CPSTIObj.Theta(j * 2)));
    end
    newCfg = tempCfgVec;
    for k = 1:Diff.num
        newCfg(Diff.pos(k)) = Diff.newVal(k);
    end
    Update.Theta = thetaShift;
    Update.tempCfgVec = newCfg;
    Diff = rmfield(Diff,'oldVal');
    Diff = rmfield(Diff,'newVal');
end
