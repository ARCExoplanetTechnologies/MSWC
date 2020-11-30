function [optics] = initialize_optics(params)

[opticalPrescription.worksheetNum opticalPrescription.worksheetTxt] = xlsread(params.opticalPrescriptionFileName);

% initialize primary plane
disp('Initializing Optical System')

for iElem = 1:params.Nelem
   
   beamrad = opticalPrescription.worksheetNum(iElem,10);
   
   tempPlaneIn = create_opticalPlane(params.Nread, 2*beamrad, params.gridCentering);
   tempPlane = create_opticalPlane(params.N, 2*beamrad, params.gridCentering);
   
   if iElem == 1
    tempNames = fieldnames(tempPlane);
    optics.elem(iElem) = tempPlane;
   end
   
   optics.elem(iElem).z = opticalPrescription.worksheetNum(iElem,4);
   optics.elem(iElem).zprop = opticalPrescription.worksheetNum(iElem,5);
   optics.elem(iElem).elemName = opticalPrescription.worksheetTxt{iElem+1,2};
   optics.elem(iElem).elemType = opticalPrescription.worksheetTxt{iElem+1,3};
   optics.elem(iElem).propType = opticalPrescription.worksheetTxt{iElem+1,6};
   optics.elem(iElem).filename = [params.dir opticalPrescription.worksheetTxt{iElem+1,8}];
   optics.elem(iElem).filenameCalib = [params.dir opticalPrescription.worksheetTxt{iElem+1,14}];
   optics.elem(iElem).R = opticalPrescription.worksheetNum(iElem,10);
   optics.elem(iElem).surfaceRMS = opticalPrescription.worksheetNum(iElem,11);
   optics.elem(iElem).reflectivityRMS = opticalPrescription.worksheetNum(iElem,12);
   optics.elem(iElem).readWriteAbb = opticalPrescription.worksheetTxt{iElem+1,13};
   optics.elem(iElem).readWriteElem = opticalPrescription.worksheetTxt{iElem+1,15};
   
   if strcmp(optics.elem(iElem).readWriteAbb, 'Write')
        optics.elem(iElem).surfaceAberrations = phase_error_flat(params.N, params.Nalpha, params.primaryAlpha, optics.elem(iElem).surfaceRMS);
        optics.elem(iElem).reflectivityAberrations = phase_error_flat(params.N, params.Nalpha, params.primaryAlpha, optics.elem(iElem).reflectivityRMS);
        fitswrite(optics.elem(iElem).surfaceAberrations, ['aberr\elem' num2str(iElem) 'surface.fits'])
        fitswrite(optics.elem(iElem).reflectivityAberrations , ['aberr\elem' num2str(iElem) 'reflectivity.fits'])
   elseif strcmp(optics.elem(iElem).readWriteAbb, 'Read')
        optics.elem(iElem).surfaceAberrations = fitsread(['aberr\elem' num2str(iElem) 'surface.fits']);
        optics.elem(iElem).reflectivityAberrations = fitsread(['aberr\elem' num2str(iElem) 'surface.fits']);
   end
   
   for j = 1:length(tempNames)
        optics.elem(iElem) = setfield(optics.elem(iElem), tempNames{j}, getfield(tempPlane, tempNames{j}));
   end
   
   if strcmp(optics.elem(iElem).elemType, 'Mirror')
       if strcmp(optics.elem(iElem).readWriteElem, 'Write')
           optics.elem(iElem).sag = zeros(params.N,params.N);
           fitswrite(optics.elem(iElem).sag, optics.elem(iElem).filename)
       elseif strcmp(optics.elem(iElem).readWriteElem, 'Read')
            sagIn = fitsread(optics.elem(iElem).filename);
            if strcmp(params.resampleMethod, 'FTdownsample')
                downsampleFactor = params.Nread / params.N;
                optics.elem(iElem).sag = FTdownsample(sagIn,downsampleFactor);
            elseif strcmp(params.resampleMethod, 'linear')
                optics.elem(iElem).sag = interp2(tempPlaneIn.xx, tempPlaneIn.yy, sagIn, optics.elem(iElem).xx, optics.elem(iElem).yy, 'linear',0);
            end 
       end
       
       if strcmp(optics.elem(iElem).elemName(1:2), 'DM')
            DMnum = str2num(optics.elem(iElem).elemName(3:end));
            optics.elem(iElem).ifArr = generate_ifarr(params.DM(DMnum).numAct,optics.elem(iElem),params.DM(DMnum).infFuncSigma);
            optics.elem(iElem).ACT = zeros(params.DM(1).numAct,params.DM(DMnum).numAct);
            optics.elem(iElem).linACT = optics.elem(iElem).ACT;
            optics.elem(iElem).flatSag = zeros(params.N,params.N);
            optics.elem(iElem).linSag = optics.elem(iElem).flatSag;
       end
       
   elseif strcmp(optics.elem(iElem).elemType, 'Binary')
        if strcmp(optics.elem(iElem).elemName, 'Lyot Stop')
            if strcmp(optics.elem(iElem).readWriteElem, 'Write')
                optics.elem(iElem).AA = (( optics.elem(iElem).rr < optics.elem(iElem).R.*params.lyotStop.FracOut) & (optics.elem(iElem).rr > optics.elem(iElem).R.*params.lyotStop.FracInn)) + 0.0;
                fitswrite(optics.elem(iElem).AA, optics.elem(iElem).filename);
            else
                optics.elem(iElem).AA = fitsread(optics.elem(iElem).filename);
            end
        else
            AAin = fitsread(optics.elem(iElem).filename) + 0.0;
            if strcmp(params.resampleMethod, 'FTdownsample')
                downsampleFactor = params.Nread / params.N;
                optics.elem(iElem).AA = FTdownsample(AAin,downsampleFactor);
            elseif strcmp(params.resampleMethod, 'linear')
                optics.elem(iElem).AA = interp2(tempPlaneIn.xx,tempPlaneIn.yy, AAin, optics.elem(iElem).xx, optics.elem(iElem).yy, 'linear',0);
            end 
        end      
   elseif strcmp(optics.elem(iElem).elemType, 'FPM')
       if strcmp(optics.elem(iElem).readWriteElem, 'Write')
            fpm = initialize_sci(params.fpm);
            fpm.M = ((fpm.rrlD < params.fpm.OutflD) & (fpm.rrlD > params.fpm.InnflD)) + 0.0;
            fpm.Mcalib = (fpm.rrlD < params.fpm.OutflD) + 0.0;     
            fitswrite(fpm.M, optics.elem(iElem).filename);
            fitswrite(fpm.M, optics.elem(iElem).filenameCalib);
       else
            fpm = initialize_sci(params.fpm);
            fpm.M = fitsread(optics.elem(iElem).filename);
            fpm.Mcalib = fitsread(optics.elem(iElem).filenameCalib);     
       end
         
       optics.elem(iElem).M = fpm.M;
       optics.elem(iElem).Mcalib = fpm.Mcalib;
       optics.elem(iElem).xxlD = fpm.xxlD;
       optics.elem(iElem).yylD = fpm.yylD;
       optics.elem(iElem).rrlD = fpm.rrlD;
       optics.elem(iElem).xx = fpm.xx;
       optics.elem(iElem).yy = fpm.yy;
       optics.elem(iElem).x = fpm.x;
       optics.elem(iElem).y = fpm.y;
       optics.elem(iElem).xlD = fpm.xlD;
       optics.elem(iElem).ylD = fpm.ylD;
       
       fitswrite(optics.elem(iElem).M, optics.elem(iElem).filename);
       elseif strcmp(optics.elem(iElem).elemType, 'Vortex')
       if strcmp(optics.elem(iElem).readWriteElem, 'Write')
            fpm = initialize_sci(params.fpm);
            %fpm.M = ((fpm.rrlD < params.fpm.OutflD) & (fpm.rrlD > params.fpm.InnflD)) + 0.0;
            fpm.M = ((fpm.rrlD < params.fpm.OutflD) + 0.0).*exp(1i*params.fpm.charge*atan2(fpm.yylD,fpm.xxlD));
            fpm.Mcalib = (fpm.rrlD < params.fpm.OutflD) + 0.0;     
            maskFilename = optics.elem(iElem).filename;
            maskFilenameRe = [maskFilename(1:(end-5)) 're.fits'];
            maskFilenameIm = [maskFilename(1:(end-5)) 'im.fits'] ;           
            fitswrite(real(fpm.M), maskFilenameRe);
            fitswrite(imag(fpm.M), maskFilenameIm);            
            fitswrite(fpm.Mcalib, optics.elem(iElem).filenameCalib);
       else
            fpm = initialize_sci(params.fpm);
            fpm.M = fitsread(optics.elem(iElem).filename);
            fpm.Mcalib = fitsread(optics.elem(iElem).filenameCalib);     
       end
         
       optics.elem(iElem).M = fpm.M;
       optics.elem(iElem).Mcalib = fpm.Mcalib;
       optics.elem(iElem).xxlD = fpm.xxlD;
       optics.elem(iElem).yylD = fpm.yylD;
       optics.elem(iElem).rrlD = fpm.rrlD;
       optics.elem(iElem).xx = fpm.xx;
       optics.elem(iElem).yy = fpm.yy;
       optics.elem(iElem).x = fpm.x;
       optics.elem(iElem).y = fpm.y;
       optics.elem(iElem).xlD = fpm.xlD;
       optics.elem(iElem).ylD = fpm.ylD;
       
       %fitswrite(optics.elem(iElem).M, optics.elem(iElem).filename);   
    elseif strcmp(optics.elem(iElem).elemType, 'FPM-Babinet')
       if strcmp(optics.elem(iElem).readWriteElem, 'Write')
            fpm = initialize_sci(params.fpm);
            fpm.M = (fpm.rrlD > params.fpm.InnflD) + 0.0;
            fpm.Mcalib = ones(size(fpm.rrlD));     
            fitswrite(fpm.M, optics.elem(iElem).filename);
            fitswrite(fpm.M, optics.elem(iElem).filenameCalib);
       else
            fpm = initialize_sci(params.fpm);
            fpm.M = fitsread(optics.elem(iElem).filename);
            fpm.Mcalib = fitsread(optics.elem(iElem).filenameCalib);     
       end
         
       optics.elem(iElem).M = fpm.M;
       optics.elem(iElem).Mcalib = fpm.Mcalib;
       optics.elem(iElem).xxlD = fpm.xxlD;
       optics.elem(iElem).yylD = fpm.yylD;
       optics.elem(iElem).rrlD = fpm.rrlD;
       optics.elem(iElem).xx = fpm.xx;
       optics.elem(iElem).yy = fpm.yy;
       optics.elem(iElem).x = fpm.x;
       optics.elem(iElem).dx = fpm.x(2) - fpm.x(1);
       optics.elem(iElem).y = fpm.y;
       optics.elem(iElem).dy = fpm.y(2) - fpm.y(1);
       optics.elem(iElem).xlD = fpm.xlD;
       optics.elem(iElem).ylD = fpm.ylD;
      
       fitswrite(optics.elem(iElem).M, optics.elem(iElem).filename);
    elseif strcmp(optics.elem(iElem).elemType, 'FPM-Babinet-IdealCMC')
       if strcmp(optics.elem(iElem).readWriteElem, 'Write')
            t = params.fpm.t;
            fpm = initialize_sci(params.fpm);
            inds = find(fpm.rrlD < params.fpm.InnflD);
            fpm.M = ones(size(fpm.rrlD));
            fpm.M(inds) = t;
            
            fpm.Mcalib = ones(size(fpm.rrlD));     
            fitswrite(fpm.M, optics.elem(iElem).filename);
            fitswrite(fpm.Mcalib, optics.elem(iElem).filenameCalib);
       else
            fpm = initialize_sci(params.fpm);
            fpm.M = fitsread(optics.elem(iElem).filename);
            fpm.Mcalib = fitsread(optics.elem(iElem).filenameCalib);     
       end
       
                
       optics.elem(iElem).M = fpm.M;
       optics.elem(iElem).Mcalib = fpm.Mcalib;
       optics.elem(iElem).xxlD = fpm.xxlD;
       optics.elem(iElem).yylD = fpm.yylD;
       optics.elem(iElem).rrlD = fpm.rrlD;
       optics.elem(iElem).xx = fpm.xx;
       optics.elem(iElem).yy = fpm.yy;
       optics.elem(iElem).x = fpm.x;
       optics.elem(iElem).y = fpm.y;
       optics.elem(iElem).dx = fpm.x(2) - fpm.x(1);
       optics.elem(iElem).dy = fpm.y(2) - fpm.y(1);
       optics.elem(iElem).xlD = fpm.xlD;
       optics.elem(iElem).ylD = fpm.ylD;
     elseif strcmp(optics.elem(iElem).elemType, 'FPM-Babinet-CsimHexCMC')
       if strcmp(optics.elem(iElem).readWriteElem, 'Write')
            t = params.fpm.t;
            fpm = initialize_sci(params.fpm);
            inds = find(fpm.rrlD < params.fpm.InnflD);
            fpm.M = ones(size(fpm.rrlD));
            fpm.M(inds) = t;
            
            fpm.Mcalib = ones(size(fpm.rrlD));     
            fitswrite(fpm.M, optics.elem(iElem).filename);
            fitswrite(fpm.M, optics.elem(iElem).filenameCalib);
       else
            fpm = initialize_sci(params.fpm);
            cmcSagFilename = optics.elem(iElem).filename;
            cmcDataHeader = cmcSagFilename(1:(end-5)); % remove file extension
            cmcSagCalibFilename = optics.elem(iElem).filenameCalib;     
       end
       
       subPixPerCell = params.fpm.subPixPerCell;
    
       optics.elem(iElem).subPixPerCell = subPixPerCell;
       optics.elem(iElem).cmcSagFilename = cmcSagFilename;
       optics.elem(iElem).cmcSagCalibFilename = cmcSagCalibFilename;
       optics.elem(iElem).maskStruct = init_mask_geometry(cmcDataHeader, cmcSagFilename, subPixPerCell);
       optics.elem(iElem).maskCalibStruct = init_mask_geometry(cmcDataHeader, cmcSagCalibFilename, subPixPerCell);
       optics.elem(iElem).xxlD = fpm.xxlD;
       optics.elem(iElem).yylD = fpm.yylD;
       optics.elem(iElem).rrlD = fpm.rrlD;
       optics.elem(iElem).xx = fpm.xx;
       optics.elem(iElem).yy = fpm.yy;
       optics.elem(iElem).x = fpm.x;
       optics.elem(iElem).y = fpm.y;
       optics.elem(iElem).dx = fpm.x(2) - fpm.x(1);
       optics.elem(iElem).dy = fpm.y(2) - fpm.y(1);
       optics.elem(iElem).xlD = fpm.xlD;
       optics.elem(iElem).ylD = fpm.ylD;   
     elseif strcmp(optics.elem(iElem).elemType, 'FPM-CMC')
       if strcmp(optics.elem(iElem).readWriteElem, 'Write')
            fpm = initialize_sci(params.fpm);
            fpm.M = (fpm.rrlD > params.fpm.InnflD) + 0.0;
            fpm.Mcalib = ones(size(fpm.rrlD));     
            fitswrite(fpm.M, optics.elem(iElem).filename);
            fitswrite(fpm.M, optics.elem(iElem).filenameCalib);
       else
            fpm = initialize_sci(params.fpm);
            %fpm.M = fitsread(optics.elem(iElem).filename);
            %fpm.Mcalib = fitsread(optics.elem(iElem).filenameCalib);     
            sagFilename = [optics.elem(iElem).filename '_sag.fits'];
            ampFilename = [optics.elem(iElem).filename '_amp.fits'];
            sagCalibFilename = [optics.elem(iElem).filenameCalib '_sag.fits'];
            ampCalibFilename = [optics.elem(iElem).filenameCalib '_amp.fits'];

            fpm.sag = fitsread(sagFilename);
            fpm.amp = fitsread(ampFilename);
            fpm.sagCalib = fitsread(sagCalibFilename);
            fpm.ampCalib = fitsread(ampCalibFilename);
       end
                
       optics.elem(iElem).sag = fpm.sag;
       optics.elem(iElem).sagCalib = fpm.sagCalib;
       optics.elem(iElem).amp = fpm.amp;
       optics.elem(iElem).ampCalib = fpm.ampCalib;
       optics.elem(iElem).xxlD = fpm.xxlD;
       optics.elem(iElem).yylD = fpm.yylD;
       optics.elem(iElem).rrlD = fpm.rrlD;
       optics.elem(iElem).xx = fpm.xx;
       optics.elem(iElem).yy = fpm.yy;
       optics.elem(iElem).x = fpm.x;
       optics.elem(iElem).y = fpm.y;
       optics.elem(iElem).xlD = fpm.xlD;
       optics.elem(iElem).ylD = fpm.ylD;
       
       %fitswrite(angle(optics.elem(iElem).M), optics.elem(iElem).filename);
    elseif strcmp(optics.elem(iElem).elemType, 'FPM-Babinet-LambdaScaling')
       if strcmp(optics.elem(iElem).readWriteElem, 'Write')
            fpm = initialize_sci(params.fpm);
            fpm.M = (fpm.rrlD > params.fpm.InnflD) + 0.0;
            fpm.Mcalib = ones(size(fpm.rrlD));     
            fitswrite(fpm.M, optics.elem(iElem).filename);
            fitswrite(fpm.M, optics.elem(iElem).filenameCalib);
       else
            fpm = initialize_sci(params.fpm);
            fpm.M = fitsread(optics.elem(iElem).filename);
            fpm.Mcalib = fitsread(optics.elem(iElem).filenameCalib);     
       end
         
       optics.elem(iElem).M = fpm.M;
       optics.elem(iElem).Mcalib = fpm.Mcalib;
       optics.elem(iElem).xxlD = fpm.xxlD;
       optics.elem(iElem).yylD = fpm.yylD;
       optics.elem(iElem).rrlD = fpm.rrlD;
       optics.elem(iElem).xx = fpm.xx;
       optics.elem(iElem).yy = fpm.yy;
       optics.elem(iElem).x = fpm.x;
       optics.elem(iElem).y = fpm.y;
       optics.elem(iElem).xlD = fpm.xlD;
       optics.elem(iElem).ylD = fpm.ylD;
       
       fitswrite(optics.elem(iElem).M, optics.elem(iElem).filename);
  elseif strcmp(optics.elem(iElem).elemType, 'FPM-Babinet-ChromaticAmplitude')
       if strcmp(optics.elem(iElem).readWriteElem, 'Write')
            fpm = initialize_sci(params.fpm);
            fpm.M = (fpm.rrlD > params.fpm.InnflD) + 0.0;
            fpm.Mcalib = ones(size(fpm.rrlD));     
            fitswrite(fpm.M, optics.elem(iElem).filename);
            fitswrite(fpm.M, optics.elem(iElem).filenameCalib);
       else
            fpm = initialize_sci(params.fpm);
            fpm.M = fitsread(optics.elem(iElem).filename);
            fpm.Mcalib = fitsread(optics.elem(iElem).filenameCalib);     
       end
        
       optics.elem(iElem).q0 = params.fpm.q0;
       optics.elem(iElem).lambda_c = params.fpm.lambda_c;
       optics.elem(iElem).B = params.fpm.B;
       optics.elem(iElem).M = fpm.M;
       optics.elem(iElem).Mcalib = fpm.Mcalib;
       optics.elem(iElem).xxlD = fpm.xxlD;
       optics.elem(iElem).yylD = fpm.yylD;
       optics.elem(iElem).rrlD = fpm.rrlD;
       optics.elem(iElem).xx = fpm.xx;
       optics.elem(iElem).yy = fpm.yy;
       optics.elem(iElem).x = fpm.x;
       optics.elem(iElem).y = fpm.y;
       optics.elem(iElem).xlD = fpm.xlD;
       optics.elem(iElem).ylD = fpm.ylD;
       
       fitswrite(optics.elem(iElem).M, optics.elem(iElem).filename);
elseif strcmp(optics.elem(iElem).elemType, 'FPM-Babinet-ChromaticAmplitude2')
       if strcmp(optics.elem(iElem).readWriteElem, 'Write')
            fpm = initialize_sci(params.fpm);
            fpm.M = (fpm.rrlD > params.fpm.InnflD) + 0.0;
            fpm.Mcalib = ones(size(fpm.rrlD));     
            fitswrite(fpm.M, optics.elem(iElem).filename);
            fitswrite(fpm.M, optics.elem(iElem).filenameCalib);
       else
            fpm = initialize_sci(params.fpm);
            fpm.M = fitsread(optics.elem(iElem).filename);
            fpm.Mcalib = fitsread(optics.elem(iElem).filenameCalib);     
       end
        
       optics.elem(iElem).q0 = params.fpm.q0;
       optics.elem(iElem).lambda_c = params.fpm.lambda_c;
       optics.elem(iElem).B = params.fpm.B;
       optics.elem(iElem).M = fpm.M;
       optics.elem(iElem).Mcalib = fpm.Mcalib;
       optics.elem(iElem).xxlD = fpm.xxlD;
       optics.elem(iElem).yylD = fpm.yylD;
       optics.elem(iElem).rrlD = fpm.rrlD;
       optics.elem(iElem).xx = fpm.xx;
       optics.elem(iElem).yy = fpm.yy;
       optics.elem(iElem).x = fpm.x;
       optics.elem(iElem).y = fpm.y;
       optics.elem(iElem).xlD = fpm.xlD;
       optics.elem(iElem).ylD = fpm.ylD;
       
       fitswrite(optics.elem(iElem).M, optics.elem(iElem).filename);
  elseif strcmp(optics.elem(iElem).elemType, 'FPM-Babinet-ChromaticAmplitude2')
       if strcmp(optics.elem(iElem).readWriteElem, 'Write')
            fpm = initialize_sci(params.fpm);
            fpm.M = (fpm.rrlD > params.fpm.InnflD) + 0.0;
            fpm.Mcalib = ones(size(fpm.rrlD));     
            fitswrite(fpm.M, optics.elem(iElem).filename);
            fitswrite(fpm.M, optics.elem(iElem).filenameCalib);
       else
            fpm = initialize_sci(params.fpm);
            fpm.M = fitsread(optics.elem(iElem).filename);
            fpm.Mcalib = fitsread(optics.elem(iElem).filenameCalib);     
       end
        
       optics.elem(iElem).q0 = params.fpm.q0;
       optics.elem(iElem).lambda_c = params.fpm.lambda_c;
       optics.elem(iElem).B = params.fpm.B;
       optics.elem(iElem).M = fpm.M;
       optics.elem(iElem).Mcalib = fpm.Mcalib;
       optics.elem(iElem).xxlD = fpm.xxlD;
       optics.elem(iElem).yylD = fpm.yylD;
       optics.elem(iElem).rrlD = fpm.rrlD;
       optics.elem(iElem).xx = fpm.xx;
       optics.elem(iElem).yy = fpm.yy;
       optics.elem(iElem).x = fpm.x;
       optics.elem(iElem).y = fpm.y;
       optics.elem(iElem).xlD = fpm.xlD;
       optics.elem(iElem).ylD = fpm.ylD;
       
       fitswrite(optics.elem(iElem).M, optics.elem(iElem).filename);     
  elseif strcmp(optics.elem(iElem).elemType, 'FPM-Babinet-ChromApod4')
       if strcmp(optics.elem(iElem).readWriteElem, 'Write')
            fpm = initialize_sci(params.fpm);
            fpm.M = (fpm.rrlD > params.fpm.InnflD) + 0.0;
            fpm.Mcalib = ones(size(fpm.rrlD));     
            fitswrite(fpm.M, optics.elem(iElem).filename);
            fitswrite(fpm.M, optics.elem(iElem).filenameCalib);
       else
            fpm = initialize_sci(params.fpm);
            fpm.M = fitsread(optics.elem(iElem).filename);
            fpm.Mcalib = fitsread(optics.elem(iElem).filenameCalib);     
       end
       
       Apodization = find_A_iterative4(100, 1000, 260, 2.6);
        
       optics.elem(iElem).q0 = params.fpm.q0;
       optics.elem(iElem).lambda_c = params.fpm.lambda_c;
       optics.elem(iElem).B = params.fpm.B;
       optics.elem(iElem).M = fpm.M;
       optics.elem(iElem).Mcalib = fpm.Mcalib;
       optics.elem(iElem).xxlD = fpm.xxlD;
       optics.elem(iElem).yylD = fpm.yylD;
       optics.elem(iElem).rrlD = fpm.rrlD;
       optics.elem(iElem).xx = fpm.xx;
       optics.elem(iElem).yy = fpm.yy;
       optics.elem(iElem).x = fpm.x;
       optics.elem(iElem).y = fpm.y;
       optics.elem(iElem).xlD = fpm.xlD;
       optics.elem(iElem).ylD = fpm.ylD;
       optics.elem(iElem).Apodization = Apodization;
       
       fitswrite(optics.elem(iElem).M, optics.elem(iElem).filename);
    elseif strcmp(optics.elem(iElem).elemType, 'FPM-Babinet-GearMask3')
       if strcmp(optics.elem(iElem).readWriteElem, 'Write')
            fpm = initialize_sci(params.fpm);
            fpm.M = (fpm.rrlD > params.fpm.InnflD) + 0.0;
            fpm.Mcalib = ones(size(fpm.rrlD));     
            fitswrite(fpm.M, optics.elem(iElem).filename);
            fitswrite(fpm.M, optics.elem(iElem).filenameCalib);
       else
            fpm = initialize_sci(params.fpm);
            fpm.M = fitsread(optics.elem(iElem).filename);
            fpm.Mcalib = fitsread(optics.elem(iElem).filenameCalib);     
       end
        
       optics.elem(iElem).q0 = params.fpm.q0;
       optics.elem(iElem).q0_mask = params.fpm.q0_mask;
       optics.elem(iElem).lambda_c = params.fpm.lambda_c;
       optics.elem(iElem).lambda_c_mask = params.fpm.lambda_c_mask;
       optics.elem(iElem).coeff1 = params.fpm.coeff1;
       optics.elem(iElem).coeff2 = params.fpm.coeff2;
       optics.elem(iElem).B = params.fpm.B;
       optics.elem(iElem).M = fpm.M;
       optics.elem(iElem).Mcalib = fpm.Mcalib;
       optics.elem(iElem).xxlD = fpm.xxlD;
       optics.elem(iElem).yylD = fpm.yylD;
       optics.elem(iElem).rrlD = fpm.rrlD;
       optics.elem(iElem).xx = fpm.xx;
       optics.elem(iElem).yy = fpm.yy;
       optics.elem(iElem).x = fpm.x;
       optics.elem(iElem).y = fpm.y;
       optics.elem(iElem).xlD = fpm.xlD;
       optics.elem(iElem).ylD = fpm.ylD;
       
       fitswrite(optics.elem(iElem).M, optics.elem(iElem).filename);
  elseif strcmp(optics.elem(iElem).elemType, 'PreFPM')
       if strcmp(optics.elem(iElem).readWriteElem, 'Write')
            prefpm = initialize_sci(params.prefpm);
            prefpm.M = ((prefpm.rrlD < params.prefpm.OutflD) & (prefpm.rrlD > params.prefpm.InnflD)) + 0.0;
            fitswrite(prefpm.M, optics.elem(iElem).filename);
       else
            prefpm = initialize_sci(params.prefpm);
            prefpm.M = fitsread(optics.elem(iElem).filename);
       end
         
       optics.elem(iElem).M = prefpm.M;
       optics.elem(iElem).xxlD = prefpm.xxlD;
       optics.elem(iElem).yylD = prefpm.yylD;
       optics.elem(iElem).rrlD = prefpm.rrlD;
       optics.elem(iElem).xx = prefpm.xx;
       optics.elem(iElem).yy = prefpm.yy;
       optics.elem(iElem).x = prefpm.x;
       optics.elem(iElem).y = prefpm.y;
       optics.elem(iElem).xlD = prefpm.xlD;
       optics.elem(iElem).ylD = prefpm.ylD;
       
       fitswrite(optics.elem(iElem).M, optics.elem(iElem).filename);    
   elseif strcmp(optics.elem(iElem).elemType, 'PIAA')
       if (strcmp(optics.elem(iElem).elemName, 'Forward') ...
               || strcmp(optics.elem(iElem).elemName, 'ForwardPhase'))
            pup = load(optics.elem(iElem).filename);
            f = pup(:,2);
       
            pupil.N = params.N;
            pupil.D = 1; 
            pupil.aapodDiam = params.PIAA.aapodDiam;
            pupil.x = ((1:pupil.N) - 1 - pupil.N/2)/(pupil.N)*pupil.D;
            pupil.y = pupil.x;
            [pupil.xx, pupil.yy] = meshgrid(pupil.x, pupil.y);
            pupil.rr = sqrt(pupil.xx.^2 + pupil.yy.^2);
            pupil.r = (0:length(f)-1)./(length(f)-1)/2/params.PIAA.pupilStop;
            pupil.ttheta = atan2(pupil.yy, pupil.xx);
       
            rt = pupil.r/max(pupil.r);
            drt = rt(2) - rt(1);
            
            Rrt_max = rt(end); % this forces the same outer edge on PIAA1 and PIAA2
            Rrt_min = params.PIAA.obstr*Rrt_max; % this forces the same obstruction on PIAA1 and PIAA2

            %Rrt = sqrt(cumtrapz(2*f.^2.*rt*drt) / trapz(2*f.^2.*rt*drt) * (Rrt_max^2 - Rrt_min^2) + Rrt_min^2);
            Rrt = sqrt(cumtrapz(2*f'.^2.*rt*drt) / trapz(2*f'.^2.*rt*drt) * (Rrt_max^2 - Rrt_min^2) + Rrt_min^2);
       
       
            PIAA.M1 = optics.elem(iElem);
            ff = interp1(pupil.r,f,pupil.rr,'spline',0).*(pupil.rr<(pupil.D/2));
            rM2_rM1_apodpre(:,1) = pupil.r;
            rM2_rM1_apodpre(:,2) = ones(size(pupil.r));
            rM2_rM1(:,1) = rt;
            rM2_rM1(:,2) = Rrt;
            rM2_apodPIAA(:,1) = pupil.r;
            rM2_apodPIAA(:,2) = f;
       
            PIAA.M1.r = pupil.r;%pupil.r*PIAA.M1.Dx/2;
            PIAA.M1.apod = ones(size(PIAA.M1.r));
            PIAA.M1.aapod = interp1(PIAA.M1.r,PIAA.M1.apod,pupil.rr,'spline').*(pupil.rr<(pupil.aapodDiam/2));

            PIAA.M2 = PIAA.M1;
            PIAA.M2.ff = ff;
            PIAA.M2.r = rM2_rM1(:,1)*optics.elem(iElem).Dx/2;
            PIAA.M2.remap = rM2_rM1(:,2)*optics.elem(iElem).Dx/2;
            PIAA.M2.A = rM2_apodPIAA(:,2);
            C = 1;
            PIAA.M2.A = PIAA.M2.A/sqrt(C);
            PIAA.M2.AA = PIAA.M2.ff;
            PIAA.M2.rremap = interp1(PIAA.M2.r, PIAA.M2.remap, PIAA.M2.rr,'spline');
            
            % measure magnification
            xi0 = -20;
            dxi = 0.05;
            xi = xi0:dxi:-xi0;
            xxi = ones(length(xi),1)*xi;
            eeta = xxi';
            tiltArr = [0 1];

            for iTilt = 1:length(tiltArr)
                tilt = tiltArr(iTilt);
                E0 = exp(2*pi*1i*tilt*pupil.xx).*(pupil.rr<=0.5);
                Epost = PIAA.M2.AA .* interp2(PIAA.M1.xx, PIAA.M1.yy, E0, PIAA.M2.rremap.*cos(PIAA.M2.ttheta), PIAA.M2.rremap.*sin(PIAA.M2.ttheta), 'spline');
    
                Epsf = zoomFFT(Epost, length(xi), dxi, xi0, length(xi), dxi, xi0);
                EpsfRef = zoomFFT(E0, length(xi), dxi, xi0, length(xi), dxi, xi0);
    
                PSF = abs(Epsf).^2;
                PSFref = abs(EpsfRef).^2;
    
                if tilt == 0
                    Peak = max(max(PSF));
                    PeakRef = max(max(PSFref));
                end
                PSF = PSF/Peak;
                PSFref = PSFref/PeakRef;
                PSF_profile = PSF(ceil(length(xi)/2),:);
                PSFref_profile = PSFref(ceil(length(xi)/2),:);

                [Y I] = max(PSF_profile);
                PSF_peak = xi(I);
                PSF_centroid = sum(sum(xxi.*PSF))/(sum(sum(PSF)));
                PSFref_centroid = sum(sum(xxi.*PSFref))/sum(sum(PSFref));
            end
            
            magnification = PSF_centroid/tilt;
            PIAA.M = magnification;

            optics.elem(iElem).PIAA = PIAA;
       elseif strcmp(optics.elem(iElem).elemName, 'Reverse')
            optics.elem(iElem).PIAA.M2inv = PIAA.M2;
            optics.elem(iElem).PIAA.M1inv = PIAA.M1;
            optics.elem(iElem).PIAA.M1inv.rremap = interp1(PIAA.M2.remap, PIAA.M2.r, PIAA.M2.rr,'spline');
                   end
   end
   
   disp(['#' num2str(iElem) '... ',  optics.elem(iElem).elemName ': ' optics.elem(iElem).filename]);
end


end

