function [tiltArr,fracEnergy] = check_iwa(system,lambda,simOptions)
      disp('Check Inner Working Angle: ')
      tiltArr = linspace(0, 4, 50);
      
      for iTilt = 1:length(tiltArr)       
            system.sources.x = tiltArr(iTilt);
            system.sources.y = 0;
            disp(['... source tilt: ' num2str(system.sources.x, '%.3f L/d')]);
       
            [results] = simResultSources(system, simOptions, lambda, system.I00, 0)
    
            Ipsf = results.Ibr;
            relativeEnergy(iTilt) = sum(sum(Ipsf));           
       end
       
       fracEnergy = relativeEnergy/max(relativeEnergy);

       % normalization
       simOptions.calibrationFlag = 1;
       system.sources.x = 0;
       system.sources.y = 0;
       noFpmResult = simResultSources(system, simOptions, lambda, system.I00, 0);
       normEnergy = sum(sum(noFpmResult.Ibr));
       fracEnergy2 = relativeEnergy/normEnergy;
       
       figure();
       plot(tiltArr, fracEnergy2);
       xlabel('Tilt, \lambda_0 / D')
       ylabel('Relative Energy')
       hold on
     
end

