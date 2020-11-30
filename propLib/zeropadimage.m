function out = zeropadimage(I, n)

% This function pads an image with zeros enough to increase the size by a
% factor of n

N = length(I)/2;
z = zeros(N*(n - 1));
zp = zeros(N*(n - 1), N*2);
out = [z zp z; zp' I zp'; z zp z];