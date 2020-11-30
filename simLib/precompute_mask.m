function [optics] = precompute_mask(optics,params)

if isfield(params,'fpm')
    if isfield(params.fpm,'usePrePropMask')
        if params.fpm.usePrePropMask == true
            % define pupil & mask
            disp('Precomputing Mask Response')
            iElem = params.fpm.elemID;
            pupil = optics.elem(iElem-1);
            pupil.gridCentering = params.gridCentering;
            pupil.E = ones(size(pupil.xx));    
            pupil.f = params.fpm.f;
            mask.babinet = params.fpm.babinet;
            mask.usePrePropMask = false;            % for first evaluation don't use prepropmask
            mask.x = optics.elem(iElem).x;
            mask.y = optics.elem(iElem).y;
            
            
            for iLambda = 1:length(params.eval.lambdaList)
                lambda = params.eval.lambdaList(iLambda);

                if strcmp(optics.elem(iElem).elemType, 'FPM-Babinet-CsimHexCMC')
                    mask.M = fliplr(compute_mask(optics.elem(iElem).maskStruct, lambda, 1));
                    maskCalib = mask;
                    maskCalib.M = fliplr(compute_mask(optics.elem(iElem).maskCalibStruct, lambda, 1));
                elseif strcmp(optics.elem(iElem).elemType, 'FPM-CMC')
                    mask.amp = optics.elem(iElem).amp;
                    mask.sag = optics.elem(iElem).sag;
                    mask.M = compute_cmc(mask.amp, mask.sag, lambda);
                    
                    maskCalib = mask;
                    maskCalib.amp = optics.elem(iElem).ampCalib;
                    maskCalib.sag = optics.elem(iElem).sagCalib;
                    maskCalib.M = compute_cmc(maskCalib.amp, maskCalib.sag, lambda);
                else
                    mask.M = optics.elem(iElem).M;
                    maskCalib = mask;
                    maskCalib.M = optics.elem(iElem).Mcalib;
                end

                disp(['@' num2str(lambda/1e-9, '%.f') ' nm'])
                [lyotOut FFT_M_hat] = pupil_to_Lyot_general(pupil, mask, lambda);
                [lyotOutCalib FFT_M_hat_calib] = pupil_to_Lyot_general(pupil, maskCalib, lambda);
                optics.elem(params.fpm.elemID).FFT_M_hat{iLambda} = FFT_M_hat;
                optics.elem(params.fpm.elemID).FFT_M_hat_calib{iLambda} = FFT_M_hat_calib;
            end
        end
    end
end

end

