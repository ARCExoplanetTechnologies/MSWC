function [Gmat] = generate_gmat(system, sourceNum, DMnum, makeNewGmat, runNum, simOptions)

lambdaCor = system.params.lambdaCor;
sci = system.sci;
regions = system.regions;
DMelemID = system.params.DM(DMnum).elemID;
Ndm = system.params.DM(DMnum).numAct;  
ifArr = system.optics.elem(DMelemID).ifArr;
I00 = system.I00;
sourceIn = system.sources(sourceNum);
simOptions.calibrationFlag = 0;
simOptions.useab = 0;
simOptions.sourceIn = sourceIn;
simOptions.verbose = 0;
simOptions.video = 0;
simOptions.videoName = 'spcVideo';
simOptions.videoFrameRate = 10;

disp('Generating G-matrix')

units
CorEle = regions.CorEle;

if simOptions.video
 v = VideoWriter(simOptions.videoName ,'mpeg-4');
 v.FrameRate = simOptions.videoFrameRate;
 open(v)
end

Einfluence = zeros(Ndm^2,length(CorEle));
EinfluenceCell = zeros([length(lambdaCor), size(Einfluence)]);



if makeNewGmat
     for iPoke = 1:(Ndm)^2,            
        DMSweepLine = zeros(1,Ndm^2);
        DMSweepLine(iPoke) = 1e-9; % poke is in nm of surface
        DMSweep = reshape(DMSweepLine,Ndm,Ndm);
        DMShapeSweep.sag = DMactToDMsurf(DMSweep,ifArr.if);
        
        % set DM poke shape as superposition of poke and nominal shape
        system.optics.elem(DMelemID).sag = system.optics.elem(DMelemID).linSag + DMShapeSweep.sag;
        
        fprintf([system.optics.elem(DMelemID).elemName ', *' num2str(sourceNum), ', ' num2str(iPoke), ': '])
        
        for iLambda = 1:length(lambdaCor)            
            lambda = lambdaCor(iLambda);
            fprintf([num2str(lambda/nm) ' nm, '])
            
            % compute un-poked propagation to obtain the poke contribution
            if iPoke == 1
                system.optics.elem(DMelemID).sag = system.optics.elem(DMelemID).linSag; % compute nominal contribution
                nomOutput = propagate_optics(system, lambda, simOptions);
                EimNom{iLambda} = nomOutput.psfE;
                system.optics.elem(DMelemID).sag = system.optics.elem(DMelemID).linSag + DMShapeSweep.sag; % add the poke back
            end
            
            sweepOutput = propagate_optics(system, lambda, simOptions);  
            sweepDelta = sweepOutput.psfE - EimNom{iLambda};
            
            if simOptions.video
                figure(103);
                subplot(1,2,1)
                imagesc(system.sci.xlD, system.sci.ylD, log10(abs(sweepDelta).^2/I00)); axis image; colorbar; caxis([-15 -8]) 
                subplot(1,2,2)
                imagesc(system.optics.elem(DMelemID).sag); axis image;
                pause(0.2);

                frame = getframe(gcf);
                writeVideo(v,frame);
            end
            
            EinfluenceCell(iLambda,iPoke,:) = sweepDelta(CorEle)/sqrt(I00);           
            SumMap(iPoke) = sum(sum(abs(sweepDelta).^2)); 
        end             
        fprintf('\n')
        
    end
    
    Gmat.EinfluenceCell = EinfluenceCell;
    Gmat.Emap = reshape(SumMap,Ndm,Ndm);
    %save(['Gmat' num2str(runNum) '.mat'],'Gmat')
    
    if simOptions.video
        close(v);
    end
else
    load(['Gmat' num2str(runNum) '.mat'])
    
end


end

