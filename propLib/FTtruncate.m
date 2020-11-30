function FTtrunc = FTtruncate(FT, N, M);

% this function truncates FT to size N x M, correctly taking care of
% the edges at the Nyquist frequency

[Nold Mold] = size(FT);

FTtrunc = [FT(1:ceil(N/2),                1:ceil(M/2))    FT(1:ceil(N/2),                 Mold - floor(M/2) + 1:Mold);
           FT(Nold - floor(N/2) + 1:Nold, 1:ceil(M/2))    FT(Nold - floor(N/2) + 1:Nold,  Mold - floor(M/2) + 1:Mold)];

if floor(N/2) == N/2 % if N is even, fix Nyquist frequency in rows
    FTtrunc(N/2+1,:) = ([FT(N/2+1,1:ceil(M/2)) FT(N/2 + 1, Mold - floor(M/2) + 1:Mold)] + FTtrunc(N/2+1,:))/2;
end

if floor(M/2) == M/2 %if M is even, fix Nyquist frequency in columns
    FTtrunc(:,M/2+1) = ([FT(1:ceil(N/2),M/2+1); FT(Nold - floor(N/2) + 1:Nold, M/2 + 1)] + FTtrunc(:,M/2+1))/2;
end

if (floor(N/2) == N/2)*(floor(M/2) == M/2) % if both are even, fix Nyquist corner
    FTtrunc(N/2+1,M/2+1) = (  FT(N/2+1,M/2+1) + FT(N/2+1,Mold - M/2 + 1)+ FT(Nold - N/2 + 1, M/2+1)+ FT(Nold - N/2 + 1, Mold - M/2 + 1)) / 4;
end