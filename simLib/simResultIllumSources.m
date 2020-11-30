%function [results] = simResultSources(optics,sci,regions,DMshape, lambda0, lambdaVec,I00,plot_flag)
function [results] = simResultIllumSources(system,simOptions,lambdaVec,I00,plot_flag)
        sources = system.sources;
        regions = system.regions;
        sci = system.sci;
               
        if plot_flag
            h = figure();
        end
        
        Ibr =  zeros(sci.N, sci.N);
        for iSource = 1:system.params.numSources
            simOptions.sourceIn = sources(iSource);
            simOptions.sourceType = 'illuminationSources';
            simOptions.iSource = iSource;
            mult = sources(iSource).mult;
            
            this_result = simResult(system,simOptions,lambdaVec,I00,plot_flag)
            Ibr = Ibr + mult*this_result.Ibr;
        end
        
        br_avg_cont = sum((Ibr(regions.score_ele)))/sum(sum(regions.ScoreMask));
        br_med_cont = median(Ibr(regions.score_ele));
        
        results.Ibr = Ibr;
        results.br_avg_cont = br_avg_cont;
        results.br_med_cont = br_med_cont;
end

