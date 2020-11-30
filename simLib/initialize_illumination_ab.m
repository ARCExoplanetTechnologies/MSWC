function [illumination] = initialize_illumination_ab(system)

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
disp('Generating Illumination')

system.illumParams = system.params;
system.opticsExport = system.optics;

if system.illumParams.illumination.downsample == true
    disp(['Setting up Illumination System through N =' num2str(system.illumParams.illumination.endElemID)])
    system.illumParams.N = system.illumParams.illumination.Nillum;
    system.illumParams.Nelem = system.illumParams.illumination.endElemID;
    system.optics = initialize_optics(system.illumParams);
end

simOptions.calibrationFlag = 1;
simOptions.useab = 1;
simOptions.sourceType = 'onaxis';
simOptions.sourceIn = system.sources(1);
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
           
	realEillumExport = interp2(system.optics.elem(simOptions.endElemID).xx, system.optics.elem(simOptions.endElemID).yy, realEillum, system.opticsExport.elem(1).xx, system.opticsExport.elem(1).yy, 'linear',0);
    imagEillumExport = interp2(system.optics.elem(simOptions.endElemID).xx, system.optics.elem(simOptions.endElemID).yy, imagEillum, system.opticsExport.elem(1).xx, system.opticsExport.elem(1).yy, 'linear',0);
    
    Eillum{iLambda}.Ein = realEillumExport + 1i*imagEillumExport;
end

illumination.lambdaList = lambdaBroadband;
illumination.Eillum = Eillum;

end
