function output = propagate_optics(system, lambda, simOptions)

sourceType = simOptions.sourceType;
useab = simOptions.useab;
thisSource = simOptions.sourceIn;
calibrationFlag = simOptions.calibrationFlag;
verbose = simOptions.verbose;
startElemID = simOptions.startElemID;
endElemID = simOptions.endElemID;

if verbose disp(['Propagating Optics: ' num2str(lambda/1e-9, '%.f') 'nm']); end;

output.elem(startElemID).Ein = 1;

if strcmp(sourceType, 'onaxis')
    output.elem(startElemID).Ein = output.elem(1).Ein;%ones(size(output.elem(1).A)).*output.elem(1).A;%.*primary.phase_error;
elseif strcmp(sourceType, 'offaxis')
    tilt = (system.params.lambdaRef / system.params.primary.D)*thisSource.x;
    tip = (system.params.lambdaRef / system.params.primary.D)*thisSource.y;
    tipTiltPhase = exp(1i*2*pi/lambda*(tilt.*system.optics.elem(1).xx + tip.*system.optics.elem(1).yy));
    output.elem(startElemID).Ein = output.elem(1).Ein.*tipTiltPhase;% = primary.E + primary.A.*primary.phase_error.*tipTiltPhase;
elseif strcmp(sourceType, 'illumination')
    iLambda = find(abs(system.illumination.lambdaList - lambda) < 1e-12);
    output.elem(startElemID).Ein = system.illumination.Eillum{iLambda}.Ein;
elseif strcmp(sourceType, 'illuminationSources')
    iLambda = find(abs(system.illumination.lambdaList - lambda) < 1e-12);
    iSource = simOptions.iSource;
    output.elem(startElemID).Ein = system.illumination.Eillum{iSource,iLambda}.Ein;
end


for iElem = startElemID:endElemID    
    if verbose disp(['#' num2str(iElem) '... ' system.optics.elem(iElem).elemName]); end
   
    if useab == 1
        output.elem(iElem).Ein = output.elem(iElem).Ein.*exp(system.optics.elem(iElem).reflectivityAberrations + 1i.*2*pi*system.optics.elem(iElem).surfaceAberrations ./ lambda);
    end
   
   if strcmp(system.optics.elem(iElem).elemType, 'Mirror')
       output.elem(iElem).Eout = output.elem(iElem).Ein.*exp(1i*4*pi/lambda*system.optics.elem(iElem).sag);       
   elseif strcmp(system.optics.elem(iElem).elemType, 'Lens')
       output.elem(iElem).Eout = output.elem(iElem).Ein.*exp(1i*2*pi/lambda*system.optics.elem(iElem).sag);
   elseif strcmp(system.optics.elem(iElem).elemType, 'Binary')
       output.elem(iElem).Eout = output.elem(iElem).Ein.*system.optics.elem(iElem).AA;
   elseif strcmp(system.optics.elem(iElem).elemType, 'FPM')  
       mask.flD = system.params.fpm.flD;%system.params.beam.f*lambda/system.optics.elem(i-1).Dx;
       mask.gridsize = system.params.fpm.FOVflD*mask.flD;
       mask.dx = mask.gridsize / system.optics.elem(iElem).N;
       mask.dy = mask.dx;
       mask.x = system.optics.elem(iElem).x;
       mask.y = system.optics.elem(iElem).y;
       mask.xx = system.optics.elem(iElem).xx;
       mask.yy = system.optics.elem(iElem).yy;
       
       mask.M = system.optics.elem(iElem).M;

       if calibrationFlag
           mask.M = system.optics.elem(iElem).Mcalib;
       end
    elseif strcmp(system.optics.elem(iElem).elemType, 'Vortex')  
       mask.flD = system.params.fpm.flD;%system.params.beam.f*lambda/system.optics.elem(i-1).Dx;
       mask.gridsize = system.params.fpm.FOVflD*mask.flD;
       mask.dx = mask.gridsize / system.optics.elem(iElem).N;
       mask.dy = mask.dx;
       mask.x = system.optics.elem(iElem).x;
       mask.y = system.optics.elem(iElem).y;
       mask.xx = system.optics.elem(iElem).xx;
       mask.yy = system.optics.elem(iElem).yy;
       
       mask.M = system.optics.elem(iElem).M;

       if calibrationFlag
           mask.M = system.optics.elem(iElem).Mcalib;
       end
     elseif strcmp(system.optics.elem(iElem).elemType, 'FPM-Babinet')  
       mask.flD = system.params.fpm.flD;%system.params.beam.f*lambda/system.optics.elem(i-1).Dx;
       mask.gridsize = system.params.fpm.FOVflD*mask.flD;
       mask.dx = mask.gridsize / system.optics.elem(iElem).N;
       mask.dy = mask.dx;
       mask.x = system.optics.elem(iElem).x;
       mask.y = system.optics.elem(iElem).y;
       mask.xx = system.optics.elem(iElem).xx;
       mask.yy = system.optics.elem(iElem).yy;
       
       mask.M = system.optics.elem(iElem).M;

       if calibrationFlag
           mask.M = system.optics.elem(iElem).Mcalib;
       end
    elseif strcmp(system.optics.elem(iElem).elemType, 'FPM-Babinet-IdealCMC')  
         maskParams.flD = system.params.fpm.f*lambda/system.params.primary.D;%system.params.beam.f*lambda/system.optics.elem(i-1).Dx;
       
         maskParams.lambdaRef = lambda;
         maskParams.samplesPerflD = system.params.fpm.samplesPerflD;
         maskParams.FOVflD = system.params.fpm.FOVflD;
       
         mask = initialize_sci(maskParams);
       
         mask.gridsize = system.params.fpm.FOVflD*mask.flD;

         %mask.M = (mask.rrlD > system.params.fpm.InnflD) + 0.0;
         inds = find(mask.rrlD < system.params.fpm.InnflD);
         mask.M = ones(size(mask.rrlD));
         mask.M(inds) = system.params.fpm.t;

       if calibrationFlag
           mask.M = system.optics.elem(iElem).Mcalib;
       end
   elseif strcmp(system.optics.elem(iElem).elemType, 'FPM-Babinet-CsimHexCMC')  
       mask.flD = system.params.fpm.flD;%system.params.beam.f*lambda/system.optics.elem(i-1).Dx;
       mask.gridsize = system.params.fpm.FOVflD*mask.flD;
       mask.dx = mask.gridsize / system.optics.elem(iElem).N;
       mask.dy = mask.dx;
       mask.x = system.optics.elem(iElem).x;
       mask.y = system.optics.elem(iElem).y;
       mask.xx = system.optics.elem(iElem).xx;
       mask.yy = system.optics.elem(iElem).yy;
       
       %mask.M = system.optics.elem(iElem).M;
       mask.M = fliplr(compute_mask(system.optics.elem(iElem).maskStruct, lambda, 1));

       if calibrationFlag
           mask.M = fliplr(compute_mask(system.optics.elem(iElem).maskCalibStruct, lambda, 1));
       end    
   elseif strcmp(system.optics.elem(iElem).elemType, 'FPM-CMC')  
       mask.flD = system.params.fpm.flD;
       mask.gridsize = system.params.fpm.FOVflD*mask.flD;
       mask.dx = mask.gridsize / system.optics.elem(iElem).N;
       mask.dy = mask.dx;
       mask.x = system.optics.elem(iElem).x;
       mask.y = system.optics.elem(iElem).y;
       mask.xx = system.optics.elem(iElem).xx;
       mask.yy = system.optics.elem(iElem).yy;
       mask.sag = system.optics.elem(iElem).sag;
       mask.sagCalib = system.optics.elem(iElem).sagCalib;
       mask.amp = system.optics.elem(iElem).amp;
       mask.ampCalib = system.optics.elem(iElem).ampCalib;

       mask.M = compute_cmc(mask.amp, mask.sag, lambda);

       if calibrationFlag
           mask.M = compute_cmc(mask.ampCalib, mask.sagCalib, lambda);
       end   
   elseif strcmp(system.optics.elem(iElem).elemType, 'FPM-Babinet-LambdaScaling')
       maskParams.flD = system.params.fpm.f*lambda/system.params.primary.D;%system.params.beam.f*lambda/system.optics.elem(i-1).Dx;
       
       maskParams.lambdaRef = lambda;
       maskParams.samplesPerflD = system.params.fpm.samplesPerflD;
       maskParams.FOVflD = system.params.fpm.FOVflD;
       
       mask = initialize_sci(maskParams);
       
       mask.gridsize = system.params.fpm.FOVflD*mask.flD;

       mask.M = (mask.rrlD > system.params.fpm.InnflD) + 0.0;
       
       if calibrationFlag
           mask.M = system.optics.elem(iElem).Mcalib;
       end
   elseif strcmp(system.optics.elem(iElem).elemType, 'FPM-Babinet-ChromApod4')
       mask.flD = system.params.fpm.flD;%system.params.beam.f*lambda/system.optics.elem(i-1).Dx;
       mask.gridsize = system.params.fpm.FOVflD*mask.flD;
       mask.dx = mask.gridsize / system.optics.elem(iElem).N;
       mask.dy = mask.dx;
       mask.x = system.optics.elem(iElem).x;
       mask.y = system.optics.elem(iElem).y;
       mask.xx = system.optics.elem(iElem).xx;
       mask.yy = system.optics.elem(iElem).yy;
       mask.rr = sqrt(system.optics.elem(iElem).xx.^2 + system.optics.elem(iElem).yy.^2);
       
       q = mask.x(ceil(end/2):end);
       q0 = system.optics.elem(iElem).q0;
       lambda_c = system.optics.elem(iElem).lambda_c;
       B = system.optics.elem(iElem).B;
       Apodization = system.optics.elem(iElem).Apodization;
       
       M1d = make_chrom_apod_mask_v4_standalone(q, q0, lambda, lambda_c, B, Apodization);
       mask.M = interp1(q,M1d,mask.rr,'pchip',1);
       
       if calibrationFlag
           mask.M = system.optics.elem(iElem).Mcalib;
       end
   elseif strcmp(system.optics.elem(iElem).elemType, 'Sci')
       exitpup = system.optics.elem(iElem-1);
       exitpup.f =  system.params.primary.f;
       exitpup.E = output.elem(iElem).Ein;
       sci = system.sci;
   end

   
   if strcmp(system.optics.elem(iElem).elemType, 'PIAA')
       if strcmp(system.optics.elem(iElem).elemName, 'Forward')
           PIAA = system.optics.elem(iElem).PIAA;
           PIAA.M1.E = output.elem(iElem).Ein.*PIAA.M1.aapod;
           PIAA.M2.E = PIAA.M2.AA .* interp2(PIAA.M1.xx, PIAA.M1.yy, PIAA.M1.E, PIAA.M2.rremap.*cos(PIAA.M2.ttheta), PIAA.M2.rremap.*sin(PIAA.M2.ttheta), 'spline');
           output.elem(iElem).Eout = PIAA.M2.E;
       elseif strcmp(system.optics.elem(iElem).elemName, 'Reverse')
           PIAA = system.optics.elem(iElem).PIAA;
           PIAA.M2inv.E = output.elem(iElem).Ein.*(PIAA.M2inv.rr < 1*PIAA.M2inv.Dx/2);
           PIAA.M1inv.E = interp2(PIAA.M2inv.xx, PIAA.M2inv.yy, 1./PIAA.M2inv.AA .* PIAA.M2inv.E, PIAA.M1inv.rremap.*cos(PIAA.M1inv.ttheta), PIAA.M1inv.rremap.*sin(PIAA.M1inv.ttheta), 'cubic');
           PIAA.M1inv.E(find(isnan(PIAA.M1inv.E))) = 0;
           PIAA.M1inv.E(find(isinf(PIAA.M1inv.E))) = 0;
           output.elem(iElem).Eout = PIAA.M1inv.E;
       elseif strcmp(system.optics.elem(iElem).elemName, 'ForwardPhase')
           PIAA = system.optics.elem(iElem).PIAA;
           PIAA.M1.E = output.elem(iElem).Ein.*PIAA.M1.aapod;
           PIAA.M2.E = PIAA.M2.AA .* interp2(PIAA.M1.xx, PIAA.M1.yy, PIAA.M1.E, PIAA.M2.rremap.*cos(PIAA.M2.ttheta), PIAA.M2.rremap.*sin(PIAA.M2.ttheta), 'spline');
           output.elem(iElem).Eout = abs(output.elem(iElem).Ein).*exp(1i*angle(PIAA.M2.E));
       end
   end
   
   if strcmp(system.optics.elem(iElem).propType, 'FresnelAS')
       output.elem(iElem+1).Ein = FresnelPropagateAS(output.elem(iElem).Eout, lambda, system.optics.elem(iElem).D/2, system.optics.elem(iElem).zprop);
   elseif strcmp(system.optics.elem(iElem).propType, 'FresnelASpad')
       padfactor = system.params.padfactor;
       output.elem(iElem+1).Ein = FresnelPropagateASpad(output.elem(iElem).Eout, lambda, system.optics.elem(iElem).D/2, system.optics.elem(iElem).zprop,padfactor);
   elseif strcmp(system.optics.elem(iElem).propType, 'PupilToLyot')       
       pupil = system.optics.elem(iElem-1);
       pupil.gridCentering = system.params.gridCentering;
       pupil.E = output.elem(iElem).Ein;    
       pupil.f = system.params.fpm.f;
       mask.lambdaList = system.params.eval.lambdaList;
       mask.babinet = system.params.fpm.babinet;
       mask.usePrePropMask = system.params.fpm.usePrePropMask;
       
       if calibrationFlag && mask.usePrePropMask
            mask.FFT_M_hat = system.optics.elem(iElem).FFT_M_hat_calib;
       elseif not(calibrationFlag) && mask.usePrePropMask
            mask.FFT_M_hat = system.optics.elem(iElem).FFT_M_hat;
       end
       
       tmp = pupil_to_Lyot_general(pupil, mask, lambda);
       output.elem(iElem).Eout = tmp.E;
       output.elem(iElem+1).Ein = tmp.E;
   elseif strcmp(system.optics.elem(iElem).propType, 'FraunhoferFocal') 
       output.psfE = Fraunhofer_focal(exitpup, sci, exitpup.f, lambda);
       output.psf = output.psfE.*conj(output.psfE);
   else
       output.elem(iElem+1).Ein = output.elem(iElem).Eout;
   end
   
   output.elem(iElem).sumIn = sum(sum(abs(output.elem(iElem).Ein).^2));
   output.lambda = lambda;
   output.simOptions = simOptions;
end

end 