%% MSWC Simulation for Roman SPC 13-dot MSWC WFOV mask
% Roman MSWC mode created by Ruslan Belikov, Eduardo Bendek, A J Riggs, and
% Dan Sirbu
% Written by Dan Sirbu, 11/27/2020
% Roman SPC WFOC Mask by A J Riggs
% Uses ZoomFFT, Chirp-Z optical propagation libraries by Rus Belikov
%
% References:
% MSWC described in:
% 
% S. Thomas, R. Belikov, and E. Bendek, "Techniques for
% High-Contrast Imaging in Multi-Star Systems I: Multi-Star Wavefront
% Control", ApJ, 810(1), 2015.
%
% D. Sirbu, S. Thomas, R. Belikov, and E. Bendek, "Techniques for
% High-Contrast Imaging in Multi-Star Systems II: Multi-Star Wavefront
% Control", ApJ, 849(2), 2017.
%
% Roman-implementation paper:
% D. Sirbu,  R. Belikov,  E. Bendek,  C. Henze,  A. J. E. Riggs,  S. Shaklan, “Multi-Star  Wavefront
% Control for the Wide-Field Infrared Survey Telescope,” in Proc. SPIE, Vol 10698, 106982F (2018)

profile on

clear all
clc
format compact
warning off all

docPref = '';
docComm = '';

libPath = 'simLib\';
utilPath = 'utilLib\';
propPath = 'propLib\';
addpath(libPath);
addpath(utilPath);
addpath(propPath);

profile on

%%  Define constants, parameters and settings
units

runNum = 1;

%% define optical system parameters, pupil, and image planes
configDir = '';
lambdaRef = 825e-9;

[runDir dmDir imagDir] = make_dirs(runNum);

disp('Initializing System Parameters')
% set config

params.dir = configDir; % optical element look-up location
params.Nelem = 8; % number of optical elements
params.Nread = 984;  % number of samples across one side of pupil planes for read-in elements
params.N = 256;
params.resampleMethod = 'FTdownsample'; % available methods are: 'FTdownsample' and 'linear'
params.lambdaRef = lambdaRef; % set reference lambda for focal plane mask
params.opticalPrescriptionFileName = 'opticalElements_romanSPC_mswc.xlsx'; % read excel file with elements
params.padfactor = 2;
params.showIntermediateFigures = true;
params.writeIntermediateFigures = true;
params.writeLogFile = true;
params.showContrastFigure = true;
params.showRelinHist = true;

% definition of primary
params.primary.D = 2*22e-3;
params.primary.fnum = 15.1; % f-number of telescope
params.primary.f = params.primary.D*params.primary.fnum; % find actual focal length
params.primarySurfaceRMS = 0; % surface errors at primary
params.primaryReflectivityRMS = 0;     % reflectivity errors at 
params.primaryAlpha = 5/2;               % aberrations ramp shape
params.flatAlphaFrac = 0.0059;         % fraction of samples FT plane for which aberrations envelope is flat
params.Nalpha = params.N*params.flatAlphaFrac;  % number of samples across FT plane for which random aberrations are flat (before flattening)

% usage of PIAA
params.usePIAA = 0;
params.useInvPIAA = 0;

% definition of PIAA
params.PIAA.pupilStop = 1;              % stopped-down PIAA fraction
params.PIAA.magnification = 1;          % assumed magnification

% beam-radius (used to size pupil planes)
params.overSizeFactor = 1;
params.beamrad = params.primary.D/2;
params.gridrad = params.overSizeFactor*params.primary.D/2;
params.gridCentering = 'cell'; %centering options are 'cell' or 'vertex'

% definition of DMs
params.numDM = 2;
params.DM(1).numAct = 48;
params.DM(1).diameter = 2*22e-3;
params.DM(1).actSpacing = params.DM(1).diameter / params.DM(1).numAct;
params.DM(1).infFuncSigma = 1.3;
params.DM(1).numActTotal = params.DM(1).numAct*params.DM(1).numAct;
params.DM(1).elemID = 3;

params.DM(2).numAct = 48;
params.DM(2).diameter = 2*22e-3;
params.DM(2).actSpacing = params.DM(1).diameter / params.DM(1).numAct;
params.DM(2).infFuncSigma = 1.3;
params.DM(2).numActTotal = params.DM(1).numAct*params.DM(1).numAct;
params.DM(2).elemID = 4;

params.numActTotal = 0;
startAct = 1;

for iDM = 1:params.numDM 
    params.DM(iDM).startAct = startAct;
    params.DM(iDM).endAct = params.DM(iDM).startAct + params.DM(iDM).numActTotal - 1;
    startAct = startAct + params.DM(iDM).numActTotal
    params.numActTotal = params.numActTotal  + params.DM(iDM).numActTotal;
    params.DM(iDM).relinNum = 1;
end

% definition of focal plane (mask) occulter
params.fpm.lambdaRef = lambdaRef;       % reference lambda for which the focal plane mask size is defined
params.fpm.InnflD = 5.6;                % focal plane mask inner radius in units of flD
params.fpm.OutflD = 20.4;                 % focal plane mask outer radius in units of flD
params.fpm.FOVflD = params.fpm.OutflD*2; % field of view used to sample focal plane mask across (twice the opening)
params.fpm.f = params.primary.f;        % focal length before the focal plane mask (set the same as for telescope/science plane)
params.fpm.flD = params.fpm.f*params.fpm.lambdaRef/params.primary.D; % define focal plane mask units of flD
params.fpm.samplesPerflD = ceil(params.N/params.fpm.FOVflD);
params.fpm.babinet = false;
params.fpm.usePrePropMask = true;
params.fpm.elemID = 6;


% definition of Lyot stop (annular)
params.lyotStop.FracInn = 0.26;         % lyot stop inner radius (as a fraction of total opening)
params.lyotStop.FracOut = 0.88;         % lyot stop outer radius (as a fraction of total opening)

% definition of science camera
params.sci.lambdaRef = lambdaRef;      % similar definition to fpm
params.sci.FOVflD = 48;
params.sci.samplesPerflD = 4;          % sampling in units of flD
params.sci.flD = params.primary.f*params.lambdaRef/params.primary.D; % definition of units of flD for the science camera

% stetch science focal plane coordinates if needed for PIAA
if params.usePIAA == 0
    if params.useInvPIAA == 0
        params.sci.flD = params.sci.flD*params.PIAA.magnification;
    end
end

% definition of scoring/control regions
params.regionType = 'rectangular';         % 'rectangular', 'annular', 'shiftAnnular'
params.ScoreMaskReg = [7 17 -5 5];        % inner and outer radius of high-contrast region
params.innerDZcutOff = 11;                  % not used here (would define inner/outer regions)
params.IWZScoreReg = [params.ScoreMaskReg(1) params.innerDZcutOff -4 4];
params.OWZScoreReg = [params.innerDZcutOff params.ScoreMaskReg(2) -4 4];
params.CorAngularSize = 360;                % angular extent in deg. of dark hole if annular
params.shiftReg = -200;                     
params.CorMaskReg = params.ScoreMaskReg;

% sources
params.numSources = 2;
params.sources(1).mult = 1;
params.sources(1).x = 0;
params.sources(1).y = 0;

params.sources(2).mult = 1;
params.sources(2).x = 112;
params.sources(2).y = 10;


% evaluation params
%params.eval.lambda0 = 575e-9;
%params.eval.lambdaList = [546:1:604]*1e-9;
params.eval.lambda0 = 825e-9;
params.eval.lambdaList = [790 808 825 843 860]*1e-9;
params.eval.broadbandProp = true;

% illumination settings
params.illumination.use = false;
params.illumination.startElemID = 1;
params.illumination.endElemID = 3;
params.illumination.downsample = true;
params.illumination.Nillum = 2048;
params.illumination.downsampleFactor = params.illumination.Nillum/params.N;
params.illumination.saveOrigIllum = true;
params.illumination.downsampleMethod = 'FTdownsample'; % choices are 'FTdownsample' and 'linear'

% correction params
params.useRelin = 1;                    % 0 - no relinearization, 1 - use relinearization
params.relinThresh = 0.1;
params.relinFlag = 0;
params.relinIterList = [];              % relinearization history tracker
params.lambdaCor = [790 808 825 843 860]*1e-9;
params.numLambdaCor = length(params.lambdaCor);
params.NiterCor = 20;
params.mu = 1e-5;


%% initialize system
system.params = params; clear params;
system.optics = initialize_optics(system.params);
system.sources = initialize_sources(system.params);
system.illumination = initialize_illumination_offaxis(system,0,0); % use same illumination as 2048x2048 (before downsampling)
system.illuminationAb = system.illumination;
system.optics = precompute_mask(system.optics,system.params);
system.sci = initialize_sci(system.params.sci);
system.regions = initialize_regions(system.sci,system.params);

%% wavelength settings
lambda0 = system.params.eval.lambda0;

if system.params.eval.broadbandProp == true
    lambdaBroadband = system.params.eval.lambdaList;
else
    lambdaBroadband = lambda0;
end

%% default simulation options
defaultSimOptions.useab = 0;
defaultSimOptions.calibrationFlag = 0;
defaultSimOptions.verbose = 1;
defaultSimOptions.sourceIn = system.sources(1);

if system.params.illumination.use == true
    defaultSimOptions.sourceType = 'illumination';
    defaultSimOptions.startElemID = system.params.illumination.endElemID + 1;
    defaultSimOptions.endElemID = system.params.Nelem;
else
    defaultSimOptions.sourceType = 'offaxis';
    defaultSimOptions.startElemID = 1;
    defaultSimOptions.endElemID = system.params.Nelem;
end


%% compute calibration PSF
disp(' ')
disp('PSF Calibration')
simOptions = defaultSimOptions;
simOptions.useab = 0;
simOptions.calibrationFlag = 1;
simOptions.sourceType = 'onaxis';
system.outputCalib = propagate_optics(system, lambda0, simOptions);
simOptions.sourceType = 'offaxis';

system.I00 = max(max(system.outputCalib.psfE.*conj(system.outputCalib.psfE)));

disp(' ')
disp('On-axis PSF')
simOptions = defaultSimOptions;
simOptions.calibrationFlag = 0;
simOptions.useab = 1;
simOptions.sourceType = 'onaxis';
system.output = propagate_optics(system, lambda0, simOptions);
simOptions.sourceType = 'offaxis';

figure(); imagesc(system.sci.xlD, system.sci.ylD, log10(system.outputCalib.psfE.*conj(system.outputCalib.psfE)./system.I00)); axis image;
set(gcf,'color','white')
title('Normalization Roman PSF', 'FontSize', 16)
xlabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
set(gcf,'color', 'white')
colorbar

figure(); imagesc(system.sci.xlD, system.sci.ylD, log10(system.output.psfE.*conj(system.output.psfE)./system.I00)); axis image;
set(gcf,'color','white')
title('Roman PSF', 'FontSize', 16)
xlabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
set(gcf,'color', 'white')
colorbar

simOptions = defaultSimOptions;
[resultsMono] = simResult(system, simOptions, lambda0, system.I00, 0)
[resultsBroadband] = simResult(system, simOptions, lambdaBroadband, system.I00, 0)
[resultsMonoMulti] = simResultSources(system, simOptions, lambda0, system.I00, 0)
[resultsBroadbandMulti] = simResultSources(system, simOptions, lambdaBroadband, system.I00, 0)


figure(); imagesc(system.sci.xlD, system.sci.ylD, log10(resultsBroadband.Ibr)); axis image; caxis([-10 -3]); colorbar; title('Normalized Intensity PSF')
set(gcf,'color','white')
title('Roman PSF - Broadband Wavelengths', 'FontSize', 16)
xlabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
set(gcf,'color', 'white')
colorbar
hold on;

figure(); imagesc(system.sci.xlD, system.sci.ylD, log10(resultsBroadbandMulti.Ibr)); axis image; caxis([-10 -3]); colorbar; title('Normalized Intensity PSF')
set(gcf,'color','white')
title('Roman PSF - Broadband Wavelengths', 'FontSize', 16)
xlabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
set(gcf,'color', 'white')
colorbar
hold on;

% generate gmatrix
simOptions = defaultSimOptions;
simOptions.useab = 0;
makeNewGmat = 1;
for iDM = 1:system.params.numDM
    for iSource = 1:system.params.numSources
        Gmat{iDM,iSource} = generate_gmat(system, iSource, iDM, makeNewGmat, runNum, simOptions);
    end
end

DMact = zeros(system.params.DM(1).numAct,system.params.DM(1).numAct);
DMnomShape.sag = DMactToDMsurf(DMact, system.optics.elem(system.params.DM(1).elemID).ifArr.if);

% initialize DM Flat & Correction shapes
for iDM = 1:system.params.numDM
    ACT{iDM} = zeros(system.params.DM(iDM).numAct,system.params.DM(iDM).numAct);
    COR{iDM}.sag = DMactToDMsurf(ACT{iDM}, system.optics.elem(system.params.DM(iDM).elemID).ifArr.if);
    system.optics.elem(system.params.DM(iDM).elemID).sag = system.optics.elem(system.params.DM(iDM).elemID).flatSag;
end

system.illuminationOnAx = system.illumination;
system.illumination = system.illuminationAb;
relinThreshHist = [];
simOptions.useab = 1;

% wavefront-control loop
for iLoop = 1:system.params.NiterCor
    tic
	fprintf(['\n Iteration ' num2str(iLoop) '\n']);
    
    % set current DM correction shape
    for iDM = 1:system.params.numDM
        system.optics.elem(system.params.DM(iDM).elemID).sag = system.optics.elem(system.params.DM(iDM).elemID).flatSag + COR{iDM}.sag;
        system.optics.elem(system.params.DM(iDM).elemID).ACT = ACT{iDM};
        fitswrite(ACT{iDM},[dmDir 'DM' num2str(iDM) 'i' num2str(iLoop) '.fits'])
    end

    intermediateResult = simResultSources(system,simOptions, lambdaBroadband, system.I00, 0);
    
    avg_cont(iLoop) = intermediateResult.br_avg_cont;
    med_cont(iLoop) = intermediateResult.br_med_cont;
    avg_cont_cor(iLoop) = intermediateResult.br_avg_cont;
    med_cont_cor(iLoop) = intermediateResult.br_med_cont;
    
    simOptions.calibrationFlag = 1;
    simOptions.useab = 0;
    
    system.illumination = system.illuminationOnAx;
    simOptions.sourceType = 'onaxis';
    system.outputSR = propagate_optics(system, lambda0,simOptions);
    simOptions.sourceType = 'offaxis';
    SRhist(iLoop) = max(max(system.outputSR.psf)) / system.I00;
    system.illumination = system.illuminationAb;
    simOptions.calibrationFlag = 0;
    
	fprintf('Mean Contrast: %.3e \n', avg_cont(iLoop));
    fprintf('Median Contrast: %.3e \n', med_cont(iLoop));
    fprintf('On-axis SR: %.3e \n',SRhist(iLoop));
   
    if system.params.showIntermediateFigures
        figure(101); imagesc(system.sci.xlD, system.sci.ylD, log10(intermediateResult.Ibr)); axis image; caxis([-9.5 -8]); colorbar; title('Normalized Intensity PSF')
        set(gcf,'color','white')
        title(['WFC Loop: ' num2str(iLoop) ', Mean Contrast: ' num2str(avg_cont(iLoop),'%.2e')], 'FontSize', 16)
        xlabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
        ylabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
        set(gcf,'color', 'white')
        colorbar
        hold on;
        caxis([-10 -3])
        pause(1)
    end
    
    if system.params.writeIntermediateFigures 
        fitswrite(intermediateResult.Ibr,[imagDir 'IM' num2str(iLoop) '.fits'])
    end
 
    if system.params.showContrastFigure
        figure(102);
        loglog(avg_cont); grid on;
        title('Contrast History')
    end

    if system.params.showRelinHist
        figure(103);
        plot(relinThreshHist); grid on;
        title('Relinearization Threshold History')
    end
    
    % check relinearization condition
    if system.params.useRelin
        for iDM = 1:system.params.numDM
            thisDMdeltaSag = system.optics.elem(system.params.DM(iDM).elemID).sag - system.optics.elem(system.params.DM(iDM).elemID).linSag;
            thisDMdelta = thisDMdeltaSag/system.params.lambdaRef*2*pi;
            thisDMdeltaMax = max(max(abs(thisDMdelta)));
            thisDM = 2*pi*system.optics.elem(system.params.DM(iDM).elemID).sag/system.params.lambdaRef;
            thisDMmax = max(max(abs(thisDM)));
            system.optics.elem(system.params.DM(iDM).elemID).deltaMaxRad = thisDMdeltaMax;
            system.optics.elem(system.params.DM(iDM).elemID).maxRad = thisDMmax;
            system.optics.elem(system.params.DM(iDM).elemID).meanRad = mean(abs(thisDM(:)));
            relinThreshHist = [relinThreshHist max(max(abs(thisDMdelta)))];

            if max(max(abs(thisDMdelta))) > system.params.relinThresh
                disp(['Relinearizing: ' num2str(thisDMdeltaMax)])
                system.params.relinFlag = 1;
                system.params.DM(iDM).relinNum = system.params.DM(iDM).relinNum + 1;
                system.params.relinIterList = [system.params.relinIterList iLoop];
                system.optics.elem(system.params.DM(iDM).elemID).linSag =  system.optics.elem(system.params.DM(iDM).elemID).sag;
                
                for iSource = 1:system.params.numSources
                    Gmat{iDM,iSource} = generate_gmat(system, iSource, iDM, makeNewGmat, runNum,simOptions);
                end
                
            end
        end
    end
    
    % create & invert G-matrix
     if ((iLoop == 1) || (system.params.relinFlag == 1))
        G = [];
        
        for iDM = 1:system.params.numDM
            G_thisDM = [];
            for iSource = 1:system.params.numSources
                for iLambda = 1:system.params.numLambdaCor
                    this_g = squeeze(Gmat{iDM,iSource}.EinfluenceCell(iLambda,:,:));
                    G_thisDM = [G_thisDM real(this_g)*system.regions.W imag(this_g)*system.regions.W];
                end
            end
            G = [G; G_thisDM];
        end
        
        [ax, ay] = size(G);
        reg = system.params.mu*diag(ones(system.params.numActTotal,1));

        ai = pinv([G reg]); 
        
        system.params.relinFlag = 0;
    end
    
    % perfect estimation
	Eest = [];
    for iSource = 1:system.params.numSources
        simOptions.sourceIn = system.params.sources(iSource);
        for iLambda = 1:system.params.numLambdaCor
            lambdaEst = system.params.lambdaCor(iLambda);
            estimateOutput = propagate_optics(system, lambdaEst, simOptions);
            this_Eest = estimateOutput.psfE(system.regions.CorEle) / sqrt(system.I00);
            Eest = [Eest; system.regions.W*real(this_Eest); system.regions.W*imag(this_Eest)];
        end
    end
        
    xreg = transpose([Eest; zeros(system.params.numActTotal,1)])*ai;
    x = -xreg * 1e-9;
     
    for iDM = 1:system.params.numDM
        xi{iDM} = x(system.params.DM(iDM).startAct:system.params.DM(iDM).endAct);
        thisACT = reshape(xi{iDM},system.params.DM(iDM).numAct, system.params.DM(iDM).numAct);
        ACT{iDM} = ACT{iDM} + thisACT;
        COR{iDM}.sag = DMactToDMsurf(ACT{iDM}, system.optics.elem(system.params.DM(iDM).elemID).ifArr.if);
        dmPVsag{iDM} = max(max(system.optics.elem(system.params.DM(iDM).elemID).sag)) - min(min(system.optics.elem(system.params.DM(iDM).elemID).sag));
        dmInds = find(abs(system.optics.elem(system.params.DM(iDM).elemID).sag) > 1e-9);
        dmRMSsag{iDM} = mean(system.optics.elem(system.params.DM(iDM).elemID).sag(dmInds));
    end
    
    loopTimeHist(iLoop) = toc;
    toc
    
    loopMetric.avg_cont = avg_cont;
    loopMetric.med_cont = med_cont;
    loopMetric.lambdas = lambdaBroadband;
    loopMetric.loopTimeHist = loopTimeHist;
    loopMetric.SRhist = SRhist;
    
    if system.params.writeLogFile
        resultsFilename = [dmDir 'results_' 'i' num2str(iLoop) '.txt']
        fid = fopen(resultsFilename,'w');
        fprintf(fid,'Mean Contrast: %.3e\r\n', avg_cont(iLoop));
        fprintf(fid,'Median Contrast: %.3e\r\n', med_cont(iLoop));
        fprintf(fid, 'On-axis SR: %.3f\r\n',SRhist(iLoop));
        fprintf(fid, 'Elapsed Time: %.3f\r\n', loopTimeHist(iLoop));
    end
    
	for iDM = 1:system.params.numDM
        dmPVsag{iDM} = max(max(system.optics.elem(system.params.DM(iDM).elemID).sag)) - min(min(system.optics.elem(system.params.DM(iDM).elemID).sag));
        dmInds = find(abs(system.optics.elem(system.params.DM(iDM).elemID).sag) > 1e-9);
        dmRMSsag{iDM} = mean(system.optics.elem(system.params.DM(iDM).elemID).sag(dmInds));
        
        fprintf(fid,'DM #%i, P-V Sag: %.3e, RMS Sag: %.3e\r\n', iDM, dmPVsag{iDM}, dmRMSsag{iDM});
    end
    
    fclose(fid)
    fclose('all')
    
 end

finalResult = simResult(system, simOptions, lambdaBroadband, system.I00, 0);
%finalResultFull = simResult(system, simOptions, lambdaBroadbandFull, system.I00, 0);

% figure(); semilogy(lambdaBroadbandFull/1e-9,finalResultFull.mono_avg_cont,'k--');
% hold on
% semilogy(lambdaBroadband/1e-9,finalResult.mono_avg_cont,'ro');
% title('Chromatic Response after EFC')
% set(gcf,'color','w')
% grid on;
% xlabel('Lambda, nm');
% ylabel('Contrast')

for iDM = 1:system.params.numDM
    figure(); imagesc(system.optics.elem(system.params.DM(iDM).elemID).sag); axis image;
    set(gcf,'color','w')
    colorbar
    PVsag = max(max(system.optics.elem(system.params.DM(iDM).elemID).sag)) - min(min(system.optics.elem(system.params.DM(iDM).elemID).sag));
    rmsInds = find(abs(system.optics.elem(system.params.DM(iDM).elemID).sag) > 1e-9);
	rmsSag = sqrt(mean(abs(system.optics.elem(system.params.DM(iDM).elemID).sag(rmsInds)).^2));
    title(['DM #', num2str(iDM), ', P-V: ', num2str(PVsag, '%.2e'), ', RMS: ', num2str(rmsSag, '%.2e')], 'FontSize', 16)
end

% figure(100); imagesc(system.sci.xlD, system.sci.ylD, log10(finalResultFull.Ibr)); axis image; caxis([-10 -3]); colorbar; title('Normalized Intensity PSF')
% set(gcf,'color','white')
% title(['Mean Contrast: ' num2str(finalResultFull.br_avg_cont,'%.2e')], 'FontSize', 16)
% xlabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
% ylabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
% set(gcf,'color', 'white')
% colorbar
% hold on;
% plot(system.regions.Cor_perimeter_x, system.regions.Cor_perimeter_y, 'Color','y')

figure(101); imagesc(system.sci.xlD, system.sci.ylD, log10(finalResult.Ibr)); axis image; caxis([-10 -3]); colorbar; title('Normalized Intensity PSF')
set(gcf,'color','white')
title(['Mean Contrast: ' num2str(finalResult.br_avg_cont,'%.2e')], 'FontSize', 16)
xlabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
set(gcf,'color', 'white')
colorbar
hold on;
plot(system.regions.Cor_perimeter_x, system.regions.Cor_perimeter_y, 'Color','y')

% simOptions.sourceType = 'offaxis';
% simOptions.sourceIn.x = 0;
% simOptions.sourceIn.y = 0;
% simOptions.startElemID = 1;
% onAxisResult = simResult(system, simOptions, lambdaBroadband, system.I00, 0);

% simOptions.sourceType = 'offaxis';
% simOptions.sourcesimOptionsIn.x = 0.004;
% simOptions.sourceIn.y = 0.004;
% tipTiltResult = simResult(system, simOptions, lambdaBroadbandFull, system.I00, 0);
% 
% simOptions.sourceType = 'offaxis';
% simOptions.sourcesimOptionsIn.x = 0.1;
% simOptions.sourceIn.y = 0.1;
% tipTiltResultLarge = simResult(system, simOptions, lambdaBroadbandFull, system.I00, 0);


% figure(102); imagesc(system.sci.xlD, system.sci.ylD, log10(tipTiltResult.Ibr)); axis image; caxis([-9.5 -8]); colorbar; title('Normalized Intensity PSF')
% set(gcf,'color','white')
% title(['Mean Contrast: ' num2str(tipTiltResult.br_avg_cont,'%.2e')], 'FontSize', 16)
% xlabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
% ylabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
% set(gcf,'color', 'white')
% colorbar
% hold on;
% plot(system.regions.Cor_perimeter_x, system.regions.Cor_perimeter_y, 'Color','y')

lambdaBroadband = system.params.eval.lambdaList;

simOptions = defaultSimOptions;

simOptionsOffAx = simOptions;
simOptionsOffAx.sourceIn = system.params.sources(2);

finalResult = simResult(system, simOptions, lambdaBroadband, system.I00, 0);
finalResultOffAx = simResult(system, simOptionsOffAx, lambdaBroadband, system.I00, 0);
finalResultSources = simResultSources(system, simOptions, lambdaBroadband, system.I00, 0);



for iDM = 1:system.params.numDM
    figure(); imagesc(system.optics.elem(system.params.DM(iDM).elemID).sag); axis image;
    set(gcf,'color','w')
    colorbar
    PVsag = max(max(system.optics.elem(system.params.DM(iDM).elemID).sag)) - min(min(system.optics.elem(system.params.DM(iDM).elemID).sag));
    title(['DM #', num2str(iDM), ', P-V Sag: ' num2str(PVsag, '%.2e')], 'FontSize', 16)
end

figure(101); imagesc(system.sci.xlD, system.sci.ylD, log10(finalResult.Ibr)); axis image; caxis([-9.5 -8]); colorbar; title('Normalized Intensity PSF')
set(gcf,'color','white')
title(['On-Axis, Mean Contrast: ' num2str(finalResult.br_avg_cont,'%.2e')], 'FontSize', 16)
xlabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
set(gcf,'color', 'white')
colorbar
hold on;
caxis([-10 -4])

figure(201); imagesc(system.sci.xlD, system.sci.ylD, log10(system.regions.CorMask.*finalResult.Ibr)); axis image; caxis([-9.5 -8]); colorbar; title('Normalized Intensity PSF')
set(gcf,'color','white')
title(['On-Axis, Mean Contrast: ' num2str(finalResult.br_avg_cont,'%.2e')], 'FontSize', 16)
xlabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
set(gcf,'color', 'white')
colorbar
hold on;
caxis([-10 -4])
plot(system.regions.Cor_perimeter_x,system.regions.Cor_perimeter_y,'y');

figure(102); imagesc(system.sci.xlD, system.sci.ylD, log10(finalResultOffAx.Ibr)); axis image; caxis([-9.5 -8]); colorbar; title('Normalized Intensity PSF')
set(gcf,'color','white')
title(['Off-Axis, Mean Contrast: ' num2str(finalResultOffAx.br_avg_cont,'%.2e')], 'FontSize', 16)
xlabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
set(gcf,'color', 'white')
colorbar
hold on;
caxis([-10 -4])


figure(202); imagesc(system.sci.xlD, system.sci.ylD, log10(system.regions.CorMask.*finalResultOffAx.Ibr)); axis image; caxis([-9.5 -8]); colorbar; title('Normalized Intensity PSF')
set(gcf,'color','white')
title(['Off-Axis, Mean Contrast: ' num2str(finalResultOffAx.br_avg_cont,'%.2e')], 'FontSize', 16)
xlabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
set(gcf,'color', 'white')
colorbar
hold on;
caxis([-10 -4])
plot(system.regions.Cor_perimeter_x,system.regions.Cor_perimeter_y,'y');


figure(103); imagesc(system.sci.xlD, system.sci.ylD, log10(finalResultSources.Ibr)); axis image; caxis([-9.5 -8]); colorbar; title('Normalized Intensity PSF')
set(gcf,'color','white')
title(['Multi-Star, Mean Contrast: ' num2str(finalResultSources.br_avg_cont,'%.2e')], 'FontSize', 16)
xlabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
set(gcf,'color', 'white')
colorbar
hold on;

simOptions.calibrationFlag = 1;
system.outputSR = propagate_optics(system, lambda0, simOptions);
SR = max(max(system.outputSR.psf)) / system.I00;
simOptions.calibrationFlag = 0;


simOptions.calibrationFlag = 1;
simOptions.sourceType = 'onaxis';
system.outputSR = propagate_optics(system, lambda0, simOptions);
simOptions.sourceType = 'offaxis';
SR = max(max(system.outputSR.psf)) / system.I00;
simOptions.calibrationFlag = 0;

% savedata
for iDM = 1:system.params.numDM
    saveData.ACT{iDM} = ACT{iDM};
    saveData.COR{iDM}.sag = COR{iDM}.sag;
end

saveData.finalResultCor = finalResult;
saveData.initialResult = resultsBroadband;
saveData.initialResultMono = resultsMono;
saveData.system.optics = system.optics;
saveData.avg_cont = avg_cont;
saveData.system.params = system.params;
saveData.system.sci = system.sci;
saveData.system.regions = system.regions;
saveData.system.outputCalib = system.outputCalib;
saveData.system.I00 = system.I00;
saveData.simOptions = simOptions;
saveData.relinThreshHist = relinThreshHist;
saveData.SRhist = SRhist;

saveFilename = [runDir 'saveData.mat']
save(saveFilename,'saveData')



profile viewer