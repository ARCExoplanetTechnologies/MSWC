function out = FresnelPropagate(x,y,Ein, xi, eta, z,lambda)

% Direct Fresnel integration based on taking a FT (zoomFFT) of the input
% field times a quadratic phase
%
% This function requires at least 4pi*F samples per aperture, where F = a^2/(lambda*Z)). 
% It works well for longer propagations where F is smaller. For larger F,
% consider using the AS version of this routine
%
% Note: this routine ignores the e^(i*k*z) and the quadratic phase factor
% at the output. This makes it as if there is a lens at the output. (To be
% corrected in the future.)

dxi = xi(2) - xi(1);
deta = eta(2) - eta(1);
dx = x(2) - x(1);
dy = y(2) - y(1);

Ax = exp(2*pi*i * dx * xi(1)/(lambda*z));
Wx = exp(-2*pi*i * dx * dxi/(lambda*z));
Ay = exp(2*pi*i * dy * eta(1)/(lambda*z));
Wy = exp(-2*pi*i * dy * deta/(lambda*z));

%Quadratic phase factor
k = 2*pi/lambda;
[xx yy] = meshgrid(x,y);
QPF = exp(i*k*(yy.^2 + xx.^2)/(2*z)); % this line needs to be optimized

out = czt(czt(Ein.*QPF, length(eta), Wy, Ay).',length(xi), Wx, Ax).'.* exp(-2*pi*i * (x(1)*ones(length(eta),1)*xi + y(1)*eta'*ones(1,length(xi)))/(lambda*z)) * dx*dy /(i*lambda*z);