function [Psi] = generateCPSTI(CPSObj,Cfg_vec)
    thetaArray = zeros(CPSObj.Nh * 2,1);
    for j = 1:CPSObj.Nh
        for I = 1:2
            theta = 0;
            for i = 1:CPSObj.Nv
                if Cfg_vec(i) ~= 0
                    theta = theta + CPSObj.W((j - 1) * CPSObj.Nv * 2 * 2 + (i - 1) * 2 * 2 + (I - 1) * 2 + Cfg_vec(i));
                elseif Cfg_vec(i) == 0
                    continue;
                end
            end
            theta = theta + CPSObj.b((j - 1) * 2 + I);
            thetaArray((j - 1) * 2 + I) = theta;
        end
    end
    Psi = 1;
    for j = 1:(CPSObj.Nh)
        Psi = Psi * (1 + exp(thetaArray(j * 2 - 1))+exp(thetaArray(j * 2)));
    end
    sum_of_a = 0;
    for i = 1:(CPSObj.Nv)
        if Cfg_vec(i) ~= 0
            sum_of_a = sum_of_a + CPSObj.a((i - 1) * 2 + Cfg_vec(i));
        end
    end
    Psi = Psi * exp(sum_of_a);
end
