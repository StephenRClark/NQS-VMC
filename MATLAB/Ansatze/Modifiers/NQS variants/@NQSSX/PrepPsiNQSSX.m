% --- General NQS wave function preparation function ---

function [NQSObj] = PrepPsiNQSSX(NQSObj,Cfg)
% This function initialises the NQS ansatz structure intermediate
% information (effective angle theta) given an initial configuration.
% ---------------------------------
% Format for NQS Modifier object with square-square interaction:
% - NQSSX.Nv = number of "visible" spins.
% - NQSSX.Nh = number of "hidden" spins.
% - NQSSX.Np = number of parameters in the ansatz = 2*Nv*Nh + 2*Nv + Nh.
% - NQSSX.a = (Nv x 1) vector - visible site bias.
% - NQSSX.A = (Nv x 1) vector - visible site square bias.
% - NQSSX.b = (Nh x 1) vector - hidden site bias.
% - NQSSX.B = (Nh x 1) vector - hidden site square bias.
% - NQSSX.W = (Nh x Nv) matrix - hidden-visible linear coupling terms.
% - NQSSX.X = (Nh x Nv) matrix - hidden-visible square coupling terms.
% - NQSSX.HDim = dimension of the hidden units.
% - NQSSX.HVal = (1 x HDim) vector of hidden unit values.
% - NQSSX.Theta = (Nh x 1) vector - effective linear-hidden angles.
% - NQSSX.VisVec = (Nv x 1) vector - visible occupancies.
% - NQSSX.ThetaSq = (Nv x 1) vector - effective square-hidden angles.
% - NQSSX.NsqVec = (Nv x 1) vector - squared visible occupancies.
% ---------------------------------
% Format for Update is a struct with two fields:
% Update.Theta - vector of new effective angles ThetaP.
% Update.VisVec - vector of new visible occupancies.
% Update.ThetaSq - vector of new effective angles ThetaSqP.
% Update.NsqVec - vector of new squared visible occupancies.
% ---------------------------------

[Cfg_vec] = NQSObj.FullCfg(Cfg); % Convert Cfg structure into a spin configuration vector.
NQSObj.VisVec = Cfg_vec;
NQSObj.NsqVec = Cfg_vec.^2; % Log the initial configuration's square values.
% Compute the effective angle for this configuration state.
NQSObj.Theta = NQSObj.b + (NQSObj.W*Cfg_vec);
NQSObj.ThetaSq = NQSObj.B + (NQSObj.X*NQSObj.NsqVec);
end