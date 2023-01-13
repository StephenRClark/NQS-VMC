function [NQSUObj] = nqsm2nqsu(NQSMObj,HilbertObj)
% This function generates a NQSU object from a NQSM object. This operation
% is lossless, as NQSU contains NQSM. Hilbert needed for initialisation
% purposes.

if HilbertObj.N ~= NQSMObj.Nv
    error('Incompatible Hilbert and NQSM site numbers.')
end

% Extract Graph for later usage.
GraphObj = NQSMObj.Graph; Nsl = max(GraphObj.SLInds);

% Use Hilbert to generate list of visible values.
Dim = HilbertObj.d; 
% List original parameters of NQSM:
Alpha = NQSMObj.Alpha; Nv = NQSMObj.Nv;
a_0 = NQSMObj.av; b_0 = NQSMObj.bv; % (Nsl*3 x 1) / (Alpha x 1) 
W_0 = NQSMObj.Wm; X_0 = NQSMObj.Xm; % (Alpha x Nv) / (Alpha x Nv)

% Target shape: a_u is (Nsl*(Dim-1) x 1), arranged [sl, vd], [sl, vd+1] ...
% a_0 reshaped is [ [H] [D] [M] ], want to replace [H] with bias for [S]
% and offset others.
a_0 = reshape(a_0,Nsl,3); a_0(:,2:3) = a_0(:,2:3) - a_0(:,1); a_0(:,1) = -a_0(:,1); 
a_u = a_0(:,3)*ones(1,(Dim-1)); a_u(:,1:3) = a_0; a_u = reshape(a_u.',Nsl*(Dim-1),1);

% Theta = b in NQSM describes all v = 1, Theta = b in NQSU describes all v = 0
% Theta_U(0) = Theta_M(0) = b_m + sum(W) = b_u
b_u = b_0 + sum(W_0);

% W is holon interaction, X is multiplon, need to interleave.
% U interactions must offset W_0 from b_u, and only differ for v = 2
% Target shape: w_u is (Alpha x Nv*(Dim-1)), arranged [a, v, vd], [a, v, vd+1]
% ... ,[a, v+1, vd], ...
w_u = zeros(Alpha,Nv*(Dim-1)); w_h = W_0; w_d = X_0;
for a = 1:Alpha
    w_temp = -ones((Dim-1),1)*w_h(a,:); w_temp(2,:) = w_temp(2,:) + w_d(a,:); % Currently ((Dim-1) x Nv)
    w_u(a,:) = w_temp(:).';
end

% Reshape these parameters into single column vector.
w_u = w_u.'; Params_U = [a_u(:); b_u; w_u(:)];

% Set ModParams for new NQSU. New W terms are not initially seeded.
ModParams.a = 0; ModParams.b = 0; ModParams.W = 0; ModParams.Alpha = Alpha;
ModParams.nmag = 0; ModParams.nphs = 0;

% Create NQSU object.
NQSUObj = NQSU(HilbertObj,GraphObj,ModParams,1); 
% Set ParamCap to permit all parameters.
ParamCap = ceil(max(abs(Params_U))/5)*5; NQSUObj.ParamCap = ParamCap;
% Load parameters into new NQS.
NQSUObj = NQSUObj.ParamLoad(Params_U);

end