function CPSTIObj = updateCPSTI(CPSTIObj,dP)
Nv = CPSTIObj.Nv;
GraphObj = CPSTIObj.Graph; Ng = GraphObj.N;
BondMap = GraphObj.BondMap;
Ntr = numel(BondMap);
Alpha = CPSTIObj.Alpha;

for i = 1:CPSTIObj.Np
  dP(i) = dP(i) * CPSTIObj.OptInds(i);
end

CPSTIObj.ati = CPSTIObj.ati + dP(1:2);
CPSTIObj.bti = CPSTIObj.bti + dP((1:(Alpha * 2)) + 2);
CPSTIObj.Wti = CPSTIObj.Wti + dP(2 + Alpha * 2 + (1:(4 * Nv * Alpha)));

cap = CPSTIObj.ParamCap;

% Sanity check the values of the ansatz:
CPSTIObj.ati(isinf(CPSTIObj.ati)) = 0;
CPSTIObj.ati(isnan(CPSTIObj.ati)) = 0;
ind = abs(CPSTIObj.ati)>cap;
CPSTIObj.ati(ind) = sign(CPSTIObj.ati(ind))*cap;

CPSTIObj.bti(isinf(CPSTIObj.bti)) = 0;
CPSTIObj.bti(isnan(CPSTIObj.bti)) = 0;
ind = abs(CPSTIObj.bti)>cap;
CPSTIObj.bti(ind) = sign(CPSTIObj.bti(ind))*cap;

CPSTIObj.Wti(isinf(CPSTIObj.Wti)) = 0;
CPSTIObj.Wti(isnan(CPSTIObj.Wti)) = 0;
ind = abs(CPSTIObj.Wti)>cap;
CPSTIObj.Wti(ind) = sign(CPSTIObj.Wti(ind))*cap;

for i = 1:Nv
    CPSTIObj.a(i * 2 - 1) = CPSTIObj.ati(1);
    CPSTIObj.a(i * 2) = CPSTIObj.ati(2);
end

for i = 1:Alpha
    for j = 1:Ntr
        for k = 1:2
            CPSTIObj.b(k + (j - 1) * 2 + (i - 1) * 2 * Ntr)=CPSTIObj.bti((i - 1) * 2 + k);
        end
    end
    for j = 1:numel(BondMap)
        for k = 1:Nv
            if BondMap{j}(1+mod(k-1,Ng)) ~= 0
                I = j + (i - 1) * Ntr;
                J = BondMap{j}(1+mod(k-1,Ng)) + Ng*(ceil(k/Ng)-1);
                CPSTIObj.W((I - 1) * Nv * 4 + (J - 1) * 4 + (1:4)) = CPSTIObj.Wti((i - 1) * Nv * 4 + (k - 1) * 4 + (1:4));
            end
        end
    end
end
