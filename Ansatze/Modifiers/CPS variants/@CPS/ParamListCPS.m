function [Params] = ParamListCPS(CPSObj)
Nv = CPSObj.Nv; Nh = CPSObj.Nh;

Params = zeros(CPSObj.Np,1);

Params(1:(Nv * 2)) = CPSObj.a;
Params((Nv * 2) + (1:(Nh * 2))) = CPSObj.b;
Params((Nv * 2) + (Nh * 2) + (1:(Nh * Nv * 4))) = CPSObj.W;
end