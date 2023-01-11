function [NQSBObj] = nqsc2nqsb(NQSCObj,HilbertObj)
% This function generates a NQSB object from a NQSC object. This operation
% is lossy, as NQSB is a truncation of NQSC. Hilbert needed for 
% initialisation purposes.

if HilbertObj.N ~= NQSCObj.Nv
    error('Incompatible Hilbert and NQSC site numbers.')
end
if HilbertObj.d ~= NQSCObj.HDim
    error('Incompatible Hilbert and NQSC site dimensions.')
end

% Extract Graph for later usage.
GraphObj = NQSCObj.Graph;

% List original parameters of NQSC:
Alpha = NQSCObj.Alpha; Params_C = NQSCObj.ParamList;

% Initialise new NQSB:
ModParams.a = 0; ModParams.b = 0; ModParams.W = 0; ModParams.Alpha = Alpha;
ModParams.nmag = 0; ModParams.nphs = 0;

% Create NQSB object.
NQSBObj = NQSB(HilbertObj,GraphObj,ModParams,1);
% Set ParamCap to permit all parameters.
ParamCap = round(max(abs(Params_C))/5)*5; NQSBObj.ParamCap = ParamCap;
% Load parameters.
Params_B = Params_C(1:NQSBObj.Np);
NQSBObj = NQSBObj.ParamLoad(Params_B);

end