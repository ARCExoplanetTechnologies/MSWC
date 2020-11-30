function [sources] = initialize_extendedSources(Npt,stellarRad,limbDarkeningFlag,u);
       if nargin < 3
           limbDarkeningFlag = 0;
       elseif nargin < 4
           u = 0.6; % limbDarkeningCoefficient
       end

        % central star
       rsource = rand(Npt,1);
       tsource = 2*pi*rand(Npt,1);
       xsource = stellarRad*sqrt(rsource).*cos(tsource);
       ysource = stellarRad*sqrt(rsource).*sin(tsource);
       
       disp(['Initializing Extended Source with Npt = ' num2str(Npt) ', R = ' num2str(stellarRad,'%.2e')])
       %sources = params.sources;
       for iSource = 1:Npt
           if limbDarkeningFlag == 0
               sources(iSource).mult = 1;
           else
               sources(iSource).mult = 1-u*(1-sqrt((1 - rsource(iSource).^2)/1));
           end
               
           sources(iSource).x = xsource(iSource);
           sources(iSource).y = ysource(iSource);
       end
end

