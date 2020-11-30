function write_output(output,writePath,subFolder)

    Nelem = length(output.elem);
    
    thisDir = [writePath filesep subFolder];
    
    mkdir(thisDir)
    
    for iElem = 1:Nelem
        iElem
        thisElemBasename = [thisDir filesep 'elem' num2str(iElem)];
      
       if not(isempty(output.elem(iElem).Ein))
            fitswrite(abs(output.elem(iElem).Ein), [thisElemBasename 'EinAmp.fits']);
            fitswrite(angle(output.elem(iElem).Ein), [thisElemBasename 'EinPha.fits']);
       end
       
       if not(isempty(output.elem(iElem).Eout))
            fitswrite(abs(output.elem(iElem).Eout), [thisElemBasename 'EoutAmp.fits']);        
            fitswrite(angle(output.elem(iElem).Eout), [thisElemBasename 'EoutPha.fits']);     
       end
    end
    
    fitswrite(output.psf, [thisDir filesep 'psf.fits']);        
    
    %save([thisDir filesep 'optics.mat'], 'optics', '-v7.3')
    
end

