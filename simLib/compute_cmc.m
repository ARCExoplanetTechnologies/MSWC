function m = compute_cmc(amp, sag, lambda, doublePassFlag)

if nargin < 4
    doublePassFlag = true;
end

if doublePassFlag 
    surfFactor = 4;
else
    surfFactor = 2;
end

m =  amp.*exp(1i*surfFactor*pi/lambda*sag);



