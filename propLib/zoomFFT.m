function out = zoomFFT(Fin, Nx, dx, x0, Ny, dy, y0)

% Computes the DFT of Fin on a Nx by Ny grid starting at (x0, y0) with
% increments (dx, dy), in units of f lambda/D. Fin is assumed to be an NxN array defined
% on a square with length and width D.

lamD = 2*pi/length(Fin);

out = czt(Fin, Nx, exp(-i*dx*lamD), exp(i*y0*lamD));
out = czt(out.',Ny, exp(-i*dy*lamD), exp(i*x0*lamD)).';