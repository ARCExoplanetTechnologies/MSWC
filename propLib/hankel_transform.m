% This function evaluates the hankel transform of Fin
% Fin is assumed to be evaluated at values r
% r is assumed to be a row vector, uniform array of values starting from 0
% output is assumed to be evaluated at q, also a row vector

function out = hankel_transform(Fin, r, q);

N = length(r);
% dr(1) = (r(2) - r(1))/2;
% dr(2:N-1) = (r(3:N) - r(1:N-2))/2;
% dr(N) = r(N) - r(N-1);
dr = r(2) - r(1);

out = (2*pi*besselj(0, 2*pi*q'*r)*(Fin.*r.*dr).').';