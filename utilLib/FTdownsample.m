% This function performs an FT-based downsampling of function f by a factor
% of M
% The size of f has to be an integer multiple of 2*M

function out = FTdownsample(f, M);

F = fft2(f);

[N1 N2] = size(f);
N1d = N1/M;
N2d = N2/M;

N1ds = [1:(N1d/2) (N1 - N1d/2 + 1):N1];
N2ds = [1:(N2d/2) (N2 - N2d/2 + 1):N2];

Fd = F(N2ds,N2ds);

%fixing the middle row and column
Fd(N1d/2 + 1, :) = (F(N1d/2 + 1, N2ds) + F(N1 - N1d/2 + 1, N2ds))/2;
Fd(:, N2d/2 + 1) = (F(N1ds, N2d/2 + 1) + F(N1ds, N2 - N2d/2 + 1))/2;
Fd(N1d/2 + 1, N2d/2 + 1) = (  F(N1d/2 + 1, N2d/2 + 1) + F(N1d/2 + 1, N2 - N2d/2 + 1) + F(N1 - N1d/2 + 1, N2 - N2d/2 + 1) + F(N1 - N1d/2 + 1, N2d/2 + 1) ) / 4;
                            
out = ifft2(Fd);