function m = compute_mask(maskGeomStruct, lambda, phaseSign)
% function m = compute_mask(maskGeomStruct, lambda, phaseSign)
%
% compute the complex mask from the data in maskGeomStruct interpolated 
% for the desired wavelength specified by lambda
%
% inputs:
% maskGeomStruct: specification of mask geometry returned by
%   init_mask_geometry
% lambda: the wavelentgh used for mask interpolation.  The resulting mask
%   is valid only for this wavelength
% phaseSign: = -1 or 1.  Optional argument to allow matching with sign 
%   conventions.  If you don't get the result you want, try changing the 
%   sign.  Defaults to 1.
%

if nargin < 3
    phaseSign = 1;
end

m = sum(maskGeomStruct.nSubPix.*(exp(phaseSign * 4*pi*1i*maskGeomStruct.sagVals/lambda)), 3)./(maskGeomStruct.subPixPerCell*maskGeomStruct.subPixPerCell);

