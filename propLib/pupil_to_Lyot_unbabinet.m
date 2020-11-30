% this function uses the FT convolution method to efficiently propagate
% from the pupil plane to the Lyot plane, multiplying by a mask in-between.
% Pupil plane is first zero-padded by a factor of 2 to avoid circular
% convolution artifacts

function lyot = pupil_to_Lyot_unbabinet(pupil, mask, lambda)

       
pupil_ext.E = zeropadimage(pupil.E,2);
pupil_ext.N = pupil.N*2;
pupil_ext.D = pupil.D*2;
pupil_ext.x = ((1:pupil_ext.N) - 1 - pupil_ext.N/2)/(pupil_ext.N)*pupil_ext.D;
pupil_ext.y = pupil_ext.x;

% propagation to mask and Lyot plane
mask.M_hat = zoomFFT_realunits(mask.x, mask.y, mask.M , pupil_ext.x, pupil_ext.y, pupil.f, lambda);

FFT_M_hat = fft2((fftshift(mask.M_hat)));
FFT_pupil = fft2(pupil_ext.E);
%E = ifft2(FFT_M_hat.*FFT_pupil)/(-1i*lambda*pupil.f)*(pupil.x(2) + pupil.x(1))^2;
E = ifft2(FFT_M_hat.*FFT_pupil)/(-1i*lambda*pupil.f)*(pupil.x(2) - pupil.x(1))^2;


% lyot plane
lyot = pupil;
lyot.E = E((1:pupil.N) + pupil.N/2, (1:pupil.N) + pupil.N/2);