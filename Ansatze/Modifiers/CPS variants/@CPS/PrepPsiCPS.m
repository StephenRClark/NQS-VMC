function [CPSObj] = PrepPsiCPS(CPSObj,Cfg)
    tempCfgVec = CPSObj.FullCfg(Cfg);
    tempCfgVec = [2;2;2];
    for i=1:(CPSObj.Nv * 2)
        CPSObj.a(i) = 0.01 * i;
    end
    for i=1:(CPSObj.Nh * 2)
        CPSObj.b(i) = 0.01 * i;
    end
    for i=1:(CPSObj.Nv * CPSObj.Nh * 4)
        CPSObj.W(i) = 0.01 * i;
    end
    %should have something to convert spin cfg to no.
    for j = 1:CPSObj.Nh
        for I = 1:2
            theta = 0;
            for i = 1:CPSObj.Nv
                if tempCfgVec(i) ~= 0
                    theta = theta + CPSObj.W((j - 1) * CPSObj.Nv * 2 * 2 + (i - 1) * 2 * 2 + (I - 1) * 2 + tempCfgVec(i));
                elseif tempCfgVec(i) == 0
                    continue;
                end
            end
            theta = theta + CPSObj.b((j - 1) * 2 + I);
            CPSObj.Theta((j - 1) * 2 + I) = theta;
        end
    end
end