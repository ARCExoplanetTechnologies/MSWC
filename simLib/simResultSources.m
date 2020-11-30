%function [results] = simResultSources(optics,sci,regions,DMshape, lambda0, lambdaVec,I00,plot_flag)
function [results] = simResultSources(system,simOptions,lambdaVec,I00,plot_flag)
        sources = system.sources;
        regions = system.regions;
        sci = system.sci;
               
        if plot_flag
            h = figure();
        end
        
        Ibr =  zeros(sci.N, sci.N);
        for isource = 1:system.params.numSources
            simOptions.sourceIn = sources(isource);
            simOptions.sourceType = 'offaxis';
            mult = sources(isource).mult;
            
            this_result = simResult(system,simOptions,lambdaVec,I00,plot_flag)
            Ibr = Ibr + mult*this_result.Ibr;
            sourceResult{isource} = this_result;
        end
        
        br_avg_cont = sum((Ibr(regions.score_ele)))/sum(sum(regions.ScoreMask));
        br_med_cont = median(Ibr(regions.score_ele));
        
        results.Ibr = Ibr;
        results.br_avg_cont = br_avg_cont;
        results.br_med_cont = br_med_cont;
        results.sourceResult = sourceResult;
end

