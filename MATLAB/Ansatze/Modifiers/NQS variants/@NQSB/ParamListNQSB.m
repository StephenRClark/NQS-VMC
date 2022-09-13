% --- General NQS wave function parameter listing function ---

function [Params] = ParamListNQSB(NQSObj)
% This function lists all the associated parameters of a general NQS in a
% single vector.
% ---------------------------------
% Format for NQSB Modifier:
% - NQSB.Nv = number of "visible" spins.
% - NQSB.Nh = number of "hidden" spins.
% - NQSB.Np = number of parameters in the ansatz = Alpha + Alpha*Nv + 2*Nsl.
% - NQSB.a = (Nv x 1) vector - visible site bias.
% - NQSB.av = (Nsl x 1) vector - visible bias parameters.
% - NQSB.A = (Nv x 1) vector - visible site square bias.
% - NQSB.Av = (Nsl x 1) vector - visible square bias parameters.
% - NQSB.b = (Nh x 1) vector - hidden site bias.
% - NQSB.bv = (Alpha x 1) vector - hidden bias parameters.
% - NQSB.B = (Nh x 1) vector- hidden site square bias.
% - NQSB.Bv = (Alpha x 1) vector - hidden square bias parameters.
% - NQSB.W = (Nh x Nv) matrix - hidden-visible coupling terms.
% - NQSB.Wm = (Alpha x Nv) matrix - coupling parameters.
% - NQSB.Alpha = number of unique coupling sets or "hidden unit density".
% - NQSB.HDim = dimension of the hidden units.
% - NQSB.Theta = (Nh x 1) vector - effective angles.
% - NQSB.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
% Format for dLogp vector is a vertically concatenated stack of parameter derivatives:
% - (Nsl x 1) for d/da.
% - (Nsl x 1) for d/dA.
% - (Alpha x 1) for d/db.
% - (Alpha x 1) for d/dB
% - (Alpha*Nv x 1) for d/dW.
% ---------------------------------

% Make local copies to reduce notation in code below.
Nv = NQSObj.Nv; % Number of "visible" spins.
% Extract information on translational symmetries from Graph.
GraphObj = NQSObj.Graph; SLInds = GraphObj.SLInds;
Nsl = max(SLInds); Alpha = NQSObj.Alpha;

Params = zeros(NQSObj.Np,1);

Params(1:Nsl) = NQSObj.av;
Params((1:Nsl)+Nsl) = NQSObj.Av;
Params((1:Alpha)+2*Nsl) = NQSObj.bv;
Params((1:Alpha)+2*Nsl+Alpha) = NQSObj.Bv;
Params((1:(Alpha*Nv))+2*Nsl+2*Alpha) = reshape(NQSObj.Wm.',Alpha*Nv,1);
end