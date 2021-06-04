function [CPSObj] = prepCPSTI(CPSObj,Cfg)
    [CPSObj.tempCfgVec] = CPSObj.FullCfg(Cfg);
    for j = 1:CPSObj.Nh
        for I = 1:2
            theta = 0;
            for i = 1:CPSObj.Nv
                if CPSObj.tempCfgVec(i) ~= 0
                    theta = theta + CPSObj.W((j - 1) * CPSObj.Nv * 2 * 2 + (i - 1) * 2 * 2 + (I - 1) * 2 + CPSObj.tempCfgVec(i));
                elseif CPSObj.tempCfgVec(i) == 0
                    continue;
                end
            end
            theta = theta + CPSObj.b((j - 1) * 2 + I);
            CPSObj.Theta((j - 1) * 2 + I) = theta;
        end
    end
end
