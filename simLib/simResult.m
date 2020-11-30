function [results] = simResult(system,simOptions,lambdaVec,I00,plot_flag)
        regions = system.regions;
        sci = system.sci;
        
        EmonoCor = zeros(length(lambdaVec),sci.N,sci.N);
        ImonoCor = zeros(length(lambdaVec),sci.N,sci.N);
        IbrCor = zeros(size(sci.N, sci.N));
    
        if plot_flag
            h = figure();
        end

        for iLambda = 1:length(lambdaVec)
            lambda = lambdaVec(iLambda);
            
            this_output = propagate_optics(system, lambda, simOptions);
            thisPSF = this_output.psf/I00;
            this_EmonoCor(iLambda,:,:) = this_output.psfE;
            this_ImonoCor(iLambda,:,:) = thisPSF;
            
            EmonoCor(iLambda,:,:) = EmonoCor(iLambda,:,:) + this_EmonoCor(iLambda,:,:);
            ImonoCor(iLambda,:,:) = ImonoCor(iLambda,:,:) + this_ImonoCor(iLambda,:,:);
            IbrCor = IbrCor + squeeze(ImonoCor(iLambda,:,:));

            mono_avg_cont(iLambda) = sum((thisPSF(regions.score_ele)))/sum(sum(regions.ScoreMask));
            med_cont(iLambda) = median(thisPSF(regions.score_ele));
            
            if plot_flag
                figure(201);
                imagesc(sci.xlD, sci.ylD, log10(flipud(thisPSF)));
                legendArr{iLambda} = [num2str(lambdaVec(iLambda)/1e-9) ' nm'];
                hold on
                caxis([-10 0])
                title(['Mono. PSF with ' num2str(lambda/1e-9) ' nm'], 'FontSize', 16)
                axis image;
                xlabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
                ylabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
                set(gcf, 'color', 'white')
                set(gca, 'FontSize', 14)
                colorbar
                axis image
                pause(0.1)
            end
            
        end
        
        IbrCor = IbrCor / length(lambdaVec);
        
        if plot_flag
            figure(h)
            plot(sci.xlD, log10(squeeze(IbrCor(ceil(sci.N/2),:))), 'k--')
            title('PSF Slice', 'FontSize', 16)
            xlabel('Angle $\xi$, $\lambda_0 / D$', 'Interpreter', 'LaTeX', 'FontSize', 16)
            ylabel('Contrast, Log$_{10}$', 'Interpreter', 'LaTeX', 'FontSize', 16)
            set(gcf, 'color', 'white')
            legendArr{iLambda+1} = 'Broadband';
            legend(legendArr)
        end
       
        br_avg_cont = sum((IbrCor(regions.score_ele)))/sum(sum(regions.ScoreMask));
        br_med_cont = median(IbrCor(regions.score_ele));
        
        if plot_flag
            figure(h);
            imagesc(sci.xlD, sci.ylD, log10(flipud(IbrCor)));
            hold on
            plot(regions.Score_perimeter_x, regions.Score_perimeter_y, 'y');
            caxis([-10 0])
            title(['Broadband PSF'], 'FontSize', 16)
            axis image;
            xlabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
            ylabel('Sky Angular Separation, $\lambda_0 / D$', 'Interpreter', 'latex', 'FontSize', 16)
            set(gcf, 'color', 'white')
            set(gca, 'FontSize', 14)
            colorbar
            axis image
        end
        
        results.med_cont = med_cont;
        results.mono_avg_cont =  mono_avg_cont;
        results.br_avg_cont = br_avg_cont;
        results.br_med_cont = br_med_cont;
        results.Emono = EmonoCor;
        results.Imono = ImonoCor;
        results.Ibr = IbrCor;
        
end