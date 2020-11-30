% this function uses the FT convolution method to efficiently propagate
% from the pupil plane to the Lyot plane, multiplying by a mask in-between.
% Pupil plane is first zero-padded by a factor of 2 to avoid circular
% convolution artifacts

function [lyot FFT_M_hat] = pupil_to_Lyot_general(pupil, mask, lambda)

if isfield(pupil,'gridCentering')
    if strcmp(pupil.gridCentering, 'vertex')
        centerFactor = 1;
    elseif strcmp(pupil.gridCentering, 'cell')
        centerFactor = 1/2;
    else
        centerFactor = 1; % default to vertex-centered
    end
else
    centerFactor = 1; % default to vertex-centered
end
        
pupil_ext.E = zeropadimage(pupil.E,2);
pupil_ext.N = pupil.N*2;
pupil_ext.D = pupil.D*2;
pupil_ext.x = ((1:pupil_ext.N) - centerFactor - pupil_ext.N/2)/(pupil_ext.N)*pupil_ext.D;
pupil_ext.y = pupil_ext.x;


% propagation to mask and Lyot plane
if not(isfield(mask,'usePrePropMask'))
    if mask.babinet == true
        mask.M_hat = zoomFFT_realunits(mask.x, mask.y, mask.M - 1, pupil_ext.x, pupil_ext.y, pupil.f, lambda);
        FFT_M_hat = fft2((fftshift(mask.M_hat)));
    elseif mask.babinet == false
        mask.M_hat = zoomFFT_realunits(mask.x, mask.y, mask.M , pupil_ext.x, pupil_ext.y, pupil.f, lambda);
        FFT_M_hat = fft2((fftshift(mask.M_hat)));
    else
        error('Define mask.babinet option')
    end
elseif mask.usePrePropMask == false
    if mask.babinet == true
        mask.M_hat = zoomFFT_realunits(mask.x, mask.y, mask.M - 1, pupil_ext.x, pupil_ext.y, pupil.f, lambda);
        FFT_M_hat = fft2((fftshift(mask.M_hat)));
    elseif mask.babinet == false
        mask.M_hat = zoomFFT_realunits(mask.x, mask.y, mask.M , pupil_ext.x, pupil_ext.y, pupil.f, lambda);
        FFT_M_hat = fft2((fftshift(mask.M_hat)));
    else
        error('Define mask.babinet option')
    end
elseif mask.usePrePropMask == true  % load prePropagated mask response corresponding to current lambda
    iLambda = find(abs(mask.lambdaList - lambda) < 1e-12);
    FFT_M_hat = mask.FFT_M_hat{iLambda};
end
   
FFT_pupil = fft2(pupil_ext.E);

if mask.babinet == false
    E = ifft2(FFT_M_hat.*FFT_pupil)/(-1i*lambda*pupil.f)*(pupil.x(2) - pupil.x(1))^2;
elseif mask.babinet == true
    E = ifft2(FFT_M_hat.*FFT_pupil)/(-1i*lambda*pupil.f)*(pupil.x(2) - pupil.x(1))^2 + pupil_ext.E;
end

% lyot plane
lyot = pupil;
lyot.E = E((1:pupil.N) + pupil.N/2, (1:pupil.N) + pupil.N/2);