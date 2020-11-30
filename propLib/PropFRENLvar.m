function [Eout,x,dx] = PropFRENLvar(Ein, z, a, b, u, v, lambda, nx, ny, figure_num)

%This function finds the propagation by evaluating the Fresnel integral
%via Matrix Fourier Transforms instead of using FFT. 

% a, b = size of first plane (aperture) in meters
% z = distance of propagation in meters
% lambda = wavelength
% nxi, neta = number of points in the first plane (odd)
% u, v = size of plane 2 in meters
% nx, ny = number of points in the second plane (odd)
% Ein = electric field at at first plane

[neta, nxi] = size(Ein);

if mod(nxi,2)==0
    %disp('Even')
    dxi = a/nxi;
    dx = u/nx;
    deta = b/neta;
    dy = v/ny;
    xi = [-nxi/2:nxi/2-1]*dxi;
    eta = [-neta/2:neta/2-1]*deta;
    x = [-nx/2:nx/2-1]*dx;
    y = [-ny/2:ny/2-1]*dy;
    x_xi = x'*xi;
    y_eta = y'*eta;
else
    %disp('Odd')
    dxi = a/(nxi-1);
    dx = u/(nx-1);
    deta = b/(neta-1);
    dy = v/(ny-1);
    xi = [-(nxi-1)/2:(nxi-1)/2]*dxi;
    eta = [-(neta-1)/2:(neta-1)/2]*deta;
    x = [-(nx-1)/2:(nx-1)/2]*dx;
    y = [-(ny-1)/2:(ny-1)/2]*dy;
    x_xi = x'*xi;
    y_eta = y'*eta;
end

Distance=log(exp(xi.^2)'*exp(eta.^2));
PropFactor=exp(1i*pi*Distance/lambda/z);

if ne(figure_num,0)

figure(figure_num)
imagesc(xi,eta,real(PropFactor))
colorbar %this commands the colorbar to turn on.  off by default
axis square %forces plotter to make plot square
title(['Propagation Factor'])
xlabel('meters')

end

Ein=Ein.*PropFactor;

Eout=exp(-2*pi*1i*y_eta/lambda/z)*Ein*exp(-2*pi*1i*x_xi/lambda/z)'*dxi*deta;

DistanceNew=log(exp(x.^2)'*exp(y.^2));

Eout=exp(2*1i*pi*z/lambda)/(1i*lambda*z)*exp(1i*pi*DistanceNew/lambda/z).*Eout;
