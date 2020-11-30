function out = Fraunhofer_focal(in, out, f, lambda)

% Fraunhofer integral (f-f lens case, with e^(i*k*z) factor dropped)
% 
% Computation is accomplished via a zoomFFT (chirpZ).
%
% "in" and "out" are structures containing x, y, dx, dy, E, and N

Ax = exp(2*pi*i * in.dx * out.x(1)/(lambda*f));
Wx = exp(-2*pi*i * in.dx * out.dx/(lambda*f));
Ay = exp(2*pi*i * in.dy * out.y(1)/(lambda*f));
Wy = exp(-2*pi*i * in.dy * out.dy/(lambda*f));

out = czt(czt(in.E, out.N, Wy, Ay).',out.N, Wx, Ax).'.* exp(-2*pi*i * (in.x(1)*ones(out.N,1)*out.x + in.y(1)*out.y'*ones(1,out.N))/(lambda*f)) * in.dx*in.dy /(i*lambda*f);