function write_coronagraph(optics,writePath)

    Nelem = length(optics.elem);
    
    thisDir = [writePath filesep 'coronagraph']
    
    mkdir(thisDir)
    
    for iElem = 1:Nelem
        thisElemBasename = [thisDir filesep 'elem' num2str(iElem)];
      
        if strcmp(optics.elem(iElem).elemType, 'Mirror')
            fitswrite(optics.elem(iElem).sag, [thisElemBasename 'sag.fits'])
        elseif strcmp(optics.elem(iElem).elemType, 'Binary')
            fitswrite(optics.elem(iElem).AA, [thisElemBasename 'amp.fits'])
        elseif strcmp(optics.elem(iElem).elemType, 'FPM-Babinet-CsimHexCMC')
            fitswrite(optics.elem(iElem).maskStruct.hexNum, [thisElemBasename 'hexNum.fits']);
            fitswrite(optics.elem(iElem).maskStruct.nSubPix, [thisElemBasename 'nSubPix.fits']);
            fitswrite(optics.elem(iElem).maskStruct.sagVals, [thisElemBasename 'sagVals.fits']);
            fitswrite(optics.elem(iElem).maskStruct.sags, [thisElemBasename 'hexSags.fits']);
            fitswrite(optics.elem(iElem).xlD, [thisElemBasename 'xlD.fits']);
            fitswrite(optics.elem(iElem).ylD, [thisElemBasename 'ylD.fits']);
        end
        
        fitswrite(optics.elem(iElem).x, [thisElemBasename 'x.fits']);
        fitswrite(optics.elem(iElem).y, [thisElemBasename 'y.fits']);
        
        fid = fopen([thisElemBasename 'info.txt'],'w');
        fprintf(fid,'R = %6.4f\r\n',optics.elem(iElem).R);
        fprintf(fid,'zprop = %6.4f\r\n',optics.elem(iElem).zprop);
        fprintf(fid,'N = %i\r\n',optics.elem(iElem).N);
        fprintf(fid,'pixscale = %6.4e\r\n',optics.elem(iElem).dx);
        fprintf(fid,'elemName = %s\r\n',optics.elem(iElem).elemName);
        fprintf(fid,'elemType = %s\r\n',optics.elem(iElem).elemType);
        fprintf(fid,'propType = %s\r\n',optics.elem(iElem).propType);
        fclose(fid);
        fclose('all')
    end
    
    %save([thisDir filesep 'optics.mat'], 'optics', '-v7.3')
    
end

