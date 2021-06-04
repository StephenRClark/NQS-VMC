function dLogp = logDerivCPSTI(NQSObj,Cfg)

Nv = NQSObj.Nv;
GraphObj = NQSObj.Graph; BondMap = GraphObj.BondMap;

Ng = GraphObj.N;
Ntr = numel(BondMap);
Alpha = NQSObj.Alpha;

Cfg_vec = NQSObj.FullCfg(Cfg);
dLogp = zeros(NQSObj.Np,1);

if NQSObj.OptInds(1) == 1
    dLogp(1) = 0;
    for i=1:Nv
        if Cfg_vec(i) == 1
            dLogp(1) = dLogp(1) + 1;
        end
    end
end

if NQSObj.OptInds(2) == 1
    dLogp(2) = 0;
    for i=1:Nv
        if Cfg_vec(i) == 2
            dLogp(2) = dLogp(2) + 1;
        end
    end
end

dThetaArray = zeros(2 * Alpha * Nv,1);
for i=1:Alpha
    for j = 1:Nv
        for k = 1:2
            I = (i - 1) * Nv * 2 + (j - 1) * 2;
            dThetaArray(I + k) = exp(NQSObj.Theta(I + k)) / (1 + exp(NQSObj.Theta(I + 1)) + exp(NQSObj.Theta(I + 2)));
        end
    end
end

for i=1:Alpha
    for k = 1:2
        for l = 1:Ntr
            I = l + (i - 1) * Ntr;
            dLogp(2 + (i - 1) * 2 + k) = dLogp(2 + (i - 1) * 2 + k) + dThetaArray((I - 1) * 2 + k);
        end
    end

    for v = 1:Nv
        for b = 1:Ntr
            TInd = b + (i - 1) * Ntr;
            VInd = BondMap{b}(1 + mod(v - 1,Ng)) + Ng * (ceil(v / Ng) - 1);
            if VInd ~= 0
                PInd = 2 + 2 * Alpha + (v - 1) * 4 + (i - 1) * Nv * 4;
                if Cfg_vec(VInd) == 1
                    dLogp(PInd+1) = dLogp(PInd+1) + dThetaArray((TInd - 1) * 2 + 1);
                    dLogp(PInd+3) = dLogp(PInd+3) + dThetaArray(TInd * 2);
                end
                if Cfg_vec(VInd) == 2
                    dLogp(PInd+2) = dLogp(PInd+2) + dThetaArray((TInd - 1) * 2 + 1);
                    dLogp(PInd+4) = dLogp(PInd+4) + dThetaArray(TInd * 2);
                end
            end
        end
    end
end

dLogp(isnan(dLogp)) = 0;
dLogp(isinf(dLogp)) = 0;
end
