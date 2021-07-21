function [P_in] = ParamCheck(P_in,cap)
% Cleans up a set of parameters, checking if they are infinite, NaN or have
% real component over the cap.
P_in(isinf(P_in)) = 0; P_in(isnan(P_in)) = 0;
ind = abs(real(P_in))>cap; 
P_in(ind) = sign(real(P_in(ind)))*cap + 1i*imag(P_in(ind));
ind = abs(imag(P_in))>pi;
P_in(ind) = real(P_in(ind)) + 1i*(mod(imag(P_in(ind))+pi,2*pi)-pi);
end