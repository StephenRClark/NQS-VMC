function [NQSMObj] = nqsu2nqsm(NQSUObj,HilbertObj)
% This function generates a NQSM object from a NQSU object. This operation
% is lossy, as NQSM is a truncation of NQSU. Hilbert needed for 
% initialisation purposes.

if HilbertObj.N ~= NQSUObj.Nv
    error('Incompatible Hilbert and NQSU site numbers.')
end

% Extract Graph for later usage.
GraphObj = NQSUObj.Graph; 

% List original parameters of NQSM:
Alpha = NQSUObj.Alpha; Nv = NQSUObj.Nv;
a_0 = NQSUObj.av; b_0 = NQSUObj.bv; % (Nsl*(Dim-1) x 1) / (Alpha x 1) 
W_0 = NQSUObj.Wm; % (Alpha x Nv*(Dim-1))

% Target shape: a_m is (3*Nsl x 1), arranged [sl, H], [sl+1, H] ... [sl, D]
% Want to calculate [H], [D], [M] from a_0, then reshape.
a_h = -a_0(:,1); a_d = a_0(:,2) - a_0(:,1); a_m = a_0(:,3) - a_0(:,1);
a_m = [a_h; a_d; a_m];

% No change to hidden bias.
b_m = b_0;

% W is holon interaction, X is doublon, need to interleave.
% Want to calculate [H], [D] from w_0, then reshape.
w_m = zeros(Alpha,Nv); x_m = zeros(Alpha,Nv);
for a = 1:Alpha
    w_temp = reshape(W_0(a,:),(Dim-1),Nv);
    w_m(a,:) = -w_temp(1,:); x_m(a,:) = w_temp(2,:) - w_temp(1,:);
end

% Reshape these parameters into single column vector.
w_m = w_m.'; x_m = x_m.'; Params_M = [a_m; b_m; w_m(:); x_m(:)];

ModParams.a = 0; ModParams.b = 0; ModParams.W = 0; ModParams.Alpha = Alpha;
ModParams.nmag = 0; ModParams.nphs = 0;

% Create NQSM object.
NQSMObj = NQSM(HilbertObj,GraphObj,ModParams,1);
% Set ParamCap to permit all parameters.
ParamCap = round(max(abs(Params_M))/5)*5; NQSMObj.ParamCap = ParamCap;
% Load parameters.
NQSMObj = NQSMObj.ParamLoad(Params_M);
end