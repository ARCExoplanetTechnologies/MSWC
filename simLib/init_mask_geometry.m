function maskGeomStruct = init_mask_geometry(maskDataHeader, sagFilename, subPixPerCell)
% function maskGeomStruct = init_mask_geometry(maskDataHeader, sagFilename, subPixPerCell)
%
% initial complex hexagonal mask geometry data, to be input to 
% compute_mask
%
% inputs:
% maskDataHeader: string used to construct the mask data file names
% sagFilename: name of the desired sag file
% subPixPerCell: optional input that must match the construction of the mask
%   geometry.  Defaults to 64, which is likely correct.
%



if nargin < 3
    subPixPerCell = 64;
end

maskGeomStruct.hexNum = fitsread([maskDataHeader '_hexNum.fits']);
maskGeomStruct.nSubPix = fitsread([maskDataHeader '_nSubPix.fits']);
maskGeomStruct.sags = fitsread(sagFilename); 
maskGeomStruct.sagVals = zeros(size(maskGeomStruct.hexNum));
maskGeomStruct.subPixPerCell = subPixPerCell;
% set the 3D sag array based on the hex number
for r = 1:size(maskGeomStruct.hexNum, 1)
    for c = 1:size(maskGeomStruct.hexNum, 2)
        for s = 1:size(maskGeomStruct.hexNum, 3)
            if (maskGeomStruct.hexNum(r, c, s) > -1)
                maskGeomStruct.sagVals(r, c, s) = maskGeomStruct.sags(maskGeomStruct.hexNum(r, c, s) + 1); % hexNum is indexed from 0
            else
                maskGeomStruct.sagVals(r, c, s) = 0;
            end
        end
    end
end