function [NQSUObj] = nqsa2nqsu(NQSAObj,HilbertObj)
% This function generates a NQSU object from a NQSA object. This operation
% is lossless, as NQSU contains NQSA. Hilbert needed for initialisation
% purposes.

if HilbertObj.N ~= NQSAObj.Nv
    error('Incompatible Hilbert and NQSA site numbers.')
end

% Extract Graph for later usage.
GraphObj = NQSAObj.Graph; Nsl = max(GraphObj.SLInds);

% Use Hilbert to generate list of visible values.
Dim = HilbertObj.d; Vals = (1:(Dim-1)).';
if strcmp(HilbertObj.Type,'Spin')
    Vals = (Vals - (Dim-1)/2)*(1+mod(Dim-1,2));
end

% List original parameters of NQSA:
Alpha = NQSAObj.Alpha; Nv = NQSAObj.Nv;
a_0 = NQSAObj.av; A_0 = NQSAObj.Av; % (Nsl x 1) vectors
b_0 = NQSAObj.bv; w_0 = NQSAObj.Wm; % (Alpha x 1) and (Alpha x Nv)

% Target shape: a_u is (Nsl*(Dim-1) x 1), arranged [sl, vd], [sl, vd+1] ...
a_u = a_0*(Vals.') + A_0*(Vals.'.^2); % Currently (Nsl x Dim-1);
a_u = reshape(a_u.',Nsl*(Dim-1),1);
% No change to hidden bias.
b_u = b_0;
% Target shape: w_u is (Alpha x Nv*(Dim-1)), arranged [a, v, vd], [a, v, vd+1]
% ... ,[a, v+1, vd], ...
w_u = zeros(Alpha,Nv*(Dim-1)); 
for a = 1:Alpha
    w_u(a,:) = reshape(Vals*w_0(a,:),Nv*(Dim-1),1).';
end

% Reshape these parameters into single column vector.
a_u = a_u.'; w_u = w_u.'; Params_U = [a_u(:); b_u; w_u(:)];

ModParams.a = 0; ModParams.b = 0; ModParams.W = 0; ModParams.Alpha = Alpha;
ModParams.nmag = 0; ModParams.nphs = 0;

% Create NQSU object.
NQSUObj = NQSU(HilbertObj,GraphObj,ModParams,1);
% Set ParamCap to permit all parameters.
ParamCap = ceil(max(abs(Params_U))/5)*5; NQSUObj.ParamCap = ParamCap;
% Load parameters into new NQS.
NQSUObj = NQSUObj.ParamLoad(Params_U);

end