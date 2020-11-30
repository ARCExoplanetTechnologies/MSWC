function show_mask(maskGeomStruct, phaseSign, colorRange)
% function show_mask(maskGeomStruct, colorRange)
%
% draw the mask sags as a 2D image
%
% inputs:
% maskGeomStruct: specification of mask geometry returned by
%   init_mask_geometry
% phaseSign: = -1 or 1.  Optional argument to allow matching with sign 
%   conventions.  If you don't get the result you want, try changing the 
%   sign.  Defaults to 1.
% colorRange: optional desired color range.  Defaults to [-5e-7 5e-7],
%   which is the allowed range for the Luvior donut design
%

if nargin < 2 || isempty(phaseSign)
    phaseSign = 1;
end
if nargin < 3 || isempty(colorRange)
    colorRange = [-5e-7 5e-7];
end

figure;
imagesc(sum(phaseSign*maskGeomStruct.nSubPix.*maskGeomStruct.sagVals, 3)./(maskGeomStruct.subPixPerCell*maskGeomStruct.subPixPerCell), colorRange);
title('FPM mask sags');
axis equal;
colorbar;
