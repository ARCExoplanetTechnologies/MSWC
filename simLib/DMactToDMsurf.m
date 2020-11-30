function DMsurf = DMactToDMsurf(DMact, InfFunc);

% This function converts from DM actuators to DM surface using the
% influence function model. The algorithm uses FFTs. This implies circular
% convolutions, i.e. circular boundary conditions
% 
% DMact: array of DM actuators. It can be complex, in which case its imaginary part represents amplitude variations.
% InfFunc: influence function, evaluated at all points in the pupil
% plane, and assumed to be centered about pupil plane origin
%
% The output DMsurf will be the same size as InfFunc. DMsurf will be
% complex-valued, but if DMact and InfFunc are real, then the imaginary part of DMsurf
% will be machine-precision close to 0.

[Nact Mact] = size(DMact);
[Nsurf Msurf] = size(InfFunc);

% Evaluate the oversized version of FT of DMact, assuming
% DMact consists of delta functions at each actuator
DMactUps = upsample(DMact, ceil((Nsurf+1)/Nact));
DMactUps = upsample(DMactUps', ceil((Msurf+1)/Mact))';
FTDMactUps = fft2(DMactUps);

% this adjusts the phase of the FTDMactUps in order to account for the shift from corner to center of each actuator
[Nups Mups] = size(FTDMactUps);
fx = [(0:ceil(Mups/2)-1) (ceil(-Mups/2):-1)]/Msurf;
fy = [(0:ceil(Nups/2)-1) (ceil(-Nups/2):-1)]'/Nsurf;
fxs = ones(Nups,1)*fx;
fys = fy*ones(1,Mups);
FTDMactUps = FTDMactUps.*exp(2*pi*i*(fxs*(Msurf*(1-1/Mact)/2) + fys*(Nsurf*(1-1/Nact)/2)));

% crop the oversized portion to Nsurf x Msurf
FTDMactTrunc = FTtruncate(FTDMactUps, Nsurf, Msurf);

FTInfFunc = fft2(InfFunc);

FTDMsurf = FTDMactTrunc.*FTInfFunc;

DMsurf = real(ifft2(FTDMsurf));
