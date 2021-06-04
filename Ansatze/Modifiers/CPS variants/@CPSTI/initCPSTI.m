function [CPSTIObj] = initCPSTI(CPSTIObj,Params)

Nv = CPSTIObj.Nv;
Nh = CPSTIObj.Nh;
GraphObj = CPSTIObj.Graph; Ng = GraphObj.N;
Alpha = round(Nh/Ng); CPSTIObj.Alpha = Alpha;
CPSTIObj.Np = 4 * Nv * Alpha + 2 * Alpha + 2;
BondMap = GraphObj.BondMap;

Ntr = numel(BondMap);

Nh = Ntr * Alpha; CPSTIObj.Nh = Nh;

CPSTIObj.a = zeros(Nv * 2,1);
CPSTIObj.ati = zeros(2,1);
CPSTIObj.b = zeros(Nh * 2,1);
CPSTIObj.bti = zeros(Alpha * 2,1);
CPSTIObj.W = zeros(Nv * Nh * 4,1);
CPSTIObj.Wti = zeros(4 * Nv * Alpha,1);
CPSTIObj.Theta = zeros(Nh * 2,1);

if isfield(Params,'A') == 0
    Params.A = Params.a;
end
if isfield(Params,'B') == 0
    Params.B = Params.b;
end

CPSTIObj.ati(1) = Params.a * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
CPSTIObj.ati(2) = Params.A * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);

for i = 1:Alpha
    CPSTIObj.bti(i * 2 - 1) = Params.b * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
    CPSTIObj.bti(i * 2) = Params.B * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
end

for i = 1:(4 * Nv * Alpha)
    CPSTIObj.Wti(i) = Params.W * (1 - Params.nmag + 2 * Params.nmag * rand) * exp(2i * pi * Params.nphs * rand);
end

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

CPSTIObj.OptInds = ones(CPSTIObj.Np,1);