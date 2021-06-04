Nv = 3; Nh = 3;
cfg = [0,0,0];
a = zeros(Nv * 2,1);
b = zeros(Nh * 2,1);
W = zeros(Nv * Nh * 4,1);
for i=1:(Nv * 2)
    a(i) = 0.01 * i;
end
for i=1:(Nh * 2)
    b(i) = 0.01 * i;
end
for i=1:(Nv * Nh * 4)
   W(i) = 0.01 * i;
end
theta = getTheta(Nv,Nh,b,W,cfg);
amplitude = getAmplitude(Nv,Nh,theta,a,cfg);

function amplitude = getAmplitude(Nv,Nh,thetaArray,a,cfg)
    amplitude = 1;
    for j = 1:Nh
        amplitude = amplitude * (1 + exp(thetaArray(j * 2 - 1))+exp(thetaArray(j * 2)));
    end
    temp = 0;
    for i = 1:Nv
        if cfg(i) ~= 0
            temp = temp + a((i - 1) * 2 + cfg(i));
        elseif cfg(i) == 0
            continue;
        end
    end
    amplitude = amplitude * exp(temp);
end

function thetaArray = getTheta(Nv,Nh,b,W,cfg)
    thetaArray = zeros(Nh * 2,1);
    for j = 1:Nh
        for I = 1:2
            theta = 0;
            for i = 1:Nv
                if cfg(i) ~= 0
                    theta = theta + W((j - 1) * Nv * 2 * 2 + (i - 1) * 2 * 2 + (I - 1) * 2 + cfg(i));
                elseif cfg(i) == 0
                    continue;
                end
            end
            theta = theta + b((j - 1) * 2 + I);
            thetaArray((j - 1) * 2 + I) = theta;
        end
    end
end
