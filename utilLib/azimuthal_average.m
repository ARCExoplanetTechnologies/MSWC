function out = azimuthal_average(A, rr, r);
% azimuthal average
% computes the azimuthal average of A(x,y)
% evaluated on the points defined by the r vector
% rr is a 2D array representing the r values of the elements of A
% r is assumed to start at 0 and be uniformly sampled

dr = r(2) - r(1);

for i = 1:length(r)
    mask = (rr >= (r(i) - dr/2)).*(rr < (r(i) + dr/2));
    out(i) = sum(sum(A.*mask))/sum(sum(mask));
end