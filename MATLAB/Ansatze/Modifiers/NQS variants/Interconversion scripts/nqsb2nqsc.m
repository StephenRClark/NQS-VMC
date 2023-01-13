function [NQSCObj] = nqsb2nqsc(NQSBObj,HilbertObj)
% This function generates a NQSC object from a NQSB object. This operation
% is lossless, as NQSC expands on NQSB. Hilbert needed for initialisation
% purposes.

if HilbertObj.N ~= NQSBObj.Nv
    error('Incompatible Hilbert and NQSB site numbers.')
end
if HilbertObj.d ~= NQSBObj.HDim
    error('Incompatible Hilbert and NQSB site dimensions.')
end

% Extract Graph for later usage.
GraphObj = NQSBObj.Graph;

% List original parameters of NQSB:
Alpha = NQSBObj.Alpha; Params_B = NQSBObj.ParamList;

% Set ModParams for new NQSC. New W terms are not initially seeded.
ModParams.a = 0; ModParams.b = 0; ModParams.W = 0; ModParams.Alpha = Alpha;
ModParams.nmag = 0; ModParams.nphs = 0;

% Create NQSC object.
NQSCObj = NQSC(HilbertObj,GraphObj,ModParams,1);
% Set ParamCap to permit all parameters.
ParamCap = ceil(max(abs(Params_B))/5)*5; NQSCObj.ParamCap = ParamCap;
% Load parameters into new NQS.
Params_C = zeros(NQSCObj.Np,1); Params_C(1:NQSBObj.Np) = Params_B;
NQSCObj = NQSCObj.ParamLoad(Params_C);
end