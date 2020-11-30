% This function simulates the phase error of a mirror as 1/f^(alpha) noise
% N is the number of samples in each dimension
% alpha is the 1/f noise power
% mirror_figure is the peak-to-valley noise variation, in arbitrary units
% out is the mirror phase noise in same units as mirror_figure

function out = phase_error_flat(N, Nflat,alpha, mirror_figure)

aberr = rand(N);
FT = fft2(aberr);

x = ((0:1:N-1)-N/2)'*ones(1,N);
y = ones(N, 1)*((0:1:N-1)-N/2);

envelope = 1./(x.^2+y.^2).^(alpha/4);

inds = find((sqrt(x.^2 + y.^2) < Nflat));
envelope(inds) = (min(envelope(inds)));

%figure(); loglog(envelope(end/2,end/2:end)) <-- plot envelope

envelope = fftshift(envelope);%fftshift(1./(x.^2+y.^2).^(alpha/4));
%envelope = fftshift(1./(x.^2+y.^2).^(alpha/4));
envelope(1,1) = 0; % DC component undefined, set to 0

FT = FT.*envelope;
out = real(fftshift(ifft2(FT)));
%out = out/(max(max(out)) - min(min(out)))*mirror_figure;
out = out/sqrt(sum(sum(out.^2)))*mirror_figure*N;