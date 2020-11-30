function out = zoomFresnelFFT_realunits(x, y, Ein, xi, eta, f, lambda)

% Fraunhover integral (f-f lens case, with e^(i*k*z) factor dropped)
% 
% Computation is accomplished via a zoomFFT (chirpZ).
%
% Ein and the output are both assumed to be defined on a square symmetric
% about the origin
%
% a = diameter of aperture in meters
% 2*Npup = number of points in pupil plane
% u = size of image plane in meters(CCD)
% 2*Nimg+1 = number of points in image plane (pixels in CCD)

[xxi eeta] = meshgrid(xi,eta);
[xx yy] = meshgrid(x,y);

quadFactor = exp(1i*pi/lambda/f*(xx.^2 + yy.^2));
%phaseMultiplier = exp(1j*pi/lambda*f)/(1i*lambda*f);
%phaseMultiplier = exp(1i*2*pi/lambda*f)*exp(1i*pi/lambda/f*(xx.^2 + yy.^2)) / (1i*lambda*f);
phaseMultiplier = exp(1i*2*pi/lambda*f).*exp(1i*pi/lambda/f*(xxi.^2 + eeta.^2));
%phaseMultiplier = exp(1i*pi/lambda/f*(xx.^2 + yy.^2));

out = phaseMultiplier.*zoomFFT_realunits(x, y, quadFactor.*Ein, xi, eta, f, lambda);