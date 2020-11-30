function [illumination,illuminationOptics] = initialize_illumination_nutek_v2(system,x,y)

%%  Define constants, parameters, and settings
units

illumination = [];

if system.params.illumination.use == true

% wavelength settings
lambda = system.params.eval.lambda0;      % default wavelength
lambdaBroadband = system.params.eval.lambdaList;

%% define optical system parameters, pupil, and image planes
units;

wfirstConfigDir = '';

disp(' ')
disp('Generating Illumination (Off-axis)')

system.illumParams = system.params;
system.opticsExport = system.optics;

if system.illumParams.illumination.downsample == true
    disp(['Setting up Illumination System at N =' num2str(system.illumParams.Nelem)])
    system.illumParams.N = system.illumParams.illumination.Nillum;
    system.illumParams.Nelem = system.illumParams.illumination.endElemID;
    
    system.optics = initialize_optics(system.illumParams);
    system.optics.elem(3).surfaceAberrations = fitsread(system.params.illumination.pupilSurf);
    system.optics.elem(4).surfaceAberrations = fitsread(system.params.illumination.nuTekSag1);
    system.optics.elem(6).surfaceAberrations = fitsread(system.params.illumination.nuTekSag2);
    system.optics.elem(7).surfaceAberrations = fitsread(system.params.illumination.auxOptics);
end

simOptions.calibrationFlag = 1;
simOptions.useab = 1;
simOptions.sourceType = 'offaxis';
simOptions.sourceIn = system.sources(1);
simOptions.sourceIn.x = x;
simOptions.sourceIn.y = y;
simOptions.verbose = 1;
simOptions.startElemID = system.params.illumination.startElemID;
simOptions.endElemID = system.params.illumination.endElemID;

% compute broadband
for iLambda = 1:length(lambdaBroadband)
    thisLambda = lambdaBroadband(iLambda);
    system.output = propagate_optics(system, thisLambda, simOptions);

	thisLambda = lambdaBroadband(iLambda);
    %disp(['Converting, at Lambda=' num2str(thisLambda) ' from nin=' num2str(length(system.optics.elem(1).xx)) ' to nout=' num2str(length(system.opticsExport.elem(1).xx))])
    
    realEillum = real(system.output.elem(simOptions.endElemID).Eout);
    imagEillum = imag(system.output.elem(simOptions.endElemID).Eout);
    
    if strcmp(system.params.illumination.downsampleMethod, 'linear')
        realEillumExport = interp2(system.optics.elem(simOptions.endElemID).xx, system.optics.elem(simOptions.endElemID).yy, realEillum, system.opticsExport.elem(1).xx, system.opticsExport.elem(1).yy, 'linear',0);
        imagEillumExport = interp2(system.optics.elem(simOptions.endElemID).xx, system.optics.elem(simOptions.endElemID).yy, imagEillum, system.opticsExport.elem(1).xx, system.opticsExport.elem(1).yy, 'linear',0);
        Eillum{iLambda}.Ein = realEillumExport + 1i*imagEillumExport;
    elseif strcmp(system.params.illumination.downsampleMethod, 'FTdownsample')
        Eillum{iLambda}.Ein = FTdownsample(system.output.elem(simOptions.endElemID).Eout,system.params.illumination.downsampleFactor);
    end
    
    
    if system.params.illumination.saveOrigIllum
        Eillum{iLambda}.Eorig = realEillum + 1i*imagEillum;
    end
end

illumination.lambdaList = lambdaBroadband;
illumination.Eillum = Eillum;
illuminationOptics = system.optics; 

end

