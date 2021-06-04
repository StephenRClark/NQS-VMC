function [Psi] = PsiGenerateCPS(CPSObj,Basis)
%I actually used configuration here
%this is not normalized
    Psi = 1;
    for j = 1:(CPSObj.Nh)
        Psi = Psi * (1 + exp(CPSObj.Theta(j * 2 - 1))+exp(CPSObj.Theta(j * 2)));
    end
    sum_of_a = 0;
    for i = 1:(CPSObj.Nv)
        if Basis(i) ~= 0
            sum_of_a = sum_of_a + CPSObj.a((i - 1) * 2 + Basis(i));
        elseif Basis(i) == 0
            continue;
        end
    end
    Psi = Psi * exp(sum_of_a);
end
