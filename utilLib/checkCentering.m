clear all;
%close all;
%filename = 'pixel_centered_circle.fits';  % John Krist pupil
%filename = 'Aperture_6_v1.fits';  % Matt V1 A
%filename = 'aperture6_rot90_inscribed_pupil2048.fits';  % Matt V1 A
%filename = 'LUVOIR-A_4096.fits';  % Matt Final A
%filename = 'LUVOIR-B_4096.fits';  % Matt Final B
%filename = 'elem1amp.fits';        % V1 A for John Krist
%filename = 'aperture6_inscribed_pupil.fits';  % Matt V1 A
%filename  = 'PIAA_1p4ld_M1_2dsag_v2-n2048.fits';
%filename = 'lyotMaskv2-rot90-remapStrutsThickStop2048.fits';
filename = 'lyotStop-0p9-2048.fits';

aperture = fitsread(filename);
aperture = circshift(aperture,[1,1]);
aperture2 = aperture(2:end,2:end);
%aperture2 = circshift(aperture,[1,1]);

figure(); imagesc(aperture - fliplr(flipud(aperture))); caxis([-1 1]); colorbar
figure(); imagesc(aperture2 - fliplr(flipud(aperture2))); caxis([-1 1]); colorbar
