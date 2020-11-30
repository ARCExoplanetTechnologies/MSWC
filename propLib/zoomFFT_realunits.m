function out = zoomFFT_realunits(x,y,Ein, xi, eta, f,lambda)

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

dxi = xi(2) - xi(1);
deta = eta(2) - eta(1);
dx = x(2) - x(1);
dy = y(2) - y(1);

Ax = exp(2*pi*i * dx * xi(1)/(lambda*f));
Wx = exp(-2*pi*i * dx * dxi/(lambda*f));
Ay = exp(2*pi*i * dy * eta(1)/(lambda*f));
Wy = exp(-2*pi*i * dy * deta/(lambda*f));

% sorry about this line, I was trying to optimize...
%out = czt(czt(Ein, 2*Nimg+1, W, A).',2*Nimg+1, W, A).'.* exp(-2*pi*i * (xs(1)*ones(2*Nimg+1,1)*xis' + ys(1)*etas*ones(1,2*Nimg+1))/(lambda*f)) * dx*dy /(lambda*f);

out = czt(czt(Ein, length(eta), Wy, Ay).',length(xi), Wx, Ax).'.* exp(-2*pi*i * (x(1)*ones(length(eta),1)*xi + y(1)*eta'*ones(1,length(xi)))/(lambda*f)) * dx*dy /(i*lambda*f);
% Below are various tests

% A = exp(-pi*i * a * u /(2*Npup*lambda*f));
% W = exp(-pi*i * a * u /(2*Npup*Nimg*lambda*f));
 
% out = czt(czt(Ein, 2*Nimg, W, A).', 2*Nimg, W, A).'.* exp(pi* i * (Npup - 1/2) * a * u * (ones(2*Nimg,1)*[-Nimg:Nimg-1] + [-Nimg:Nimg-1]'*ones(1,2*Nimg)) /(2*Nimg*Npup*lambda*f)) * a^2 /(4*Npup^2*lambda*f);
% out = czt(czt(Ein, 2*Nimg, W, A).', 2*Nimg, W, A).'.* exp(i*(ones(2*Nimg,1)*[-Nimg:Nimg-1]));