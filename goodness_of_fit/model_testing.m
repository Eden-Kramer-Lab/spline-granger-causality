%% Model Test #1 - Durbin-Watson test for autocorrelation of model residuals

[notwhite, dw] = gof_residuals(model_spline);
if isempty(notwhite)
    fprintf('For the spline-Granger fit, all signals pass Durbin-Watson test, FDR corrected, q = %g\n',0.05);
else
    fprintf(2,'For spline-Granger fit, there are autocorrelated residuals, FDR corrected, q = %g, for signals(s): %s\n',0.05,num2str(notwhite));
end
[notwhite, dw] = gof_residuals(model_standard);

if isempty(notwhite)
    fprintf('For the standard-Granger fit, all signals pass Durbin-Watson test, FDR corrected, q = %g\n',0.05);
else
    fprintf(2,'For standard-Granger fit, there are autocorrelated residuals, FDR corrected, q = %g, for signals(s): %s\n',0.05,num2str(notwhite));
end

%% Model Test #2 - Grenander and Rosenblatt test of the Integrated Spectrum

[spline_fit] = grstat(model_true,model_spline);

if isempty(spline_fit.fails)
    fprintf('For the spline-Granger fit, all signals pass GR test, Bonferroni corrected, at significance %g\n',0.05);
else
    fprintf(2,'For spline-Granger fit, signal(s) fail GR test, Bonferroni corrected, at significance %g, for: %s\n',0.05,num2str(spline_fit.fails));
end


[standard_fit] = grstat(model_true,model_standard);

if isempty(standard_fit.fails)
    fprintf('For the standard-Granger fit, all signals pass GR test, Bonferroni corrected, at significance %g\n',0.05);
else
    fprintf(2,'For standard-Granger fit, signal(s) fail GR test, Bonferroni corrected, at significance %g, for: %s\n',0.05,num2str(standard_fit.fails));
end

nelectrodes = size(model_spline.network,1);

if model_true.show_all_plots
    figure;
    for i = 1:nelectrodes
        subplot(sqrt(nelectrodes),sqrt(nelectrodes),i)
        
        spline_axis = spline_fit.xaxis(i,:);
        standard_axis = standard_fit.xaxis(i,:);
        
        plot(spline_axis,spline_fit.estimate(i,:),'r','LineWidth',1.5)
        hold on
        plot(standard_axis,standard_fit.estimate(i,:),'g','LineWidth',1.5)
        plot(standard_axis,standard_fit.true(i,:),'k','LineWidth',1.5)
        
        plot(spline_axis,spline_fit.bound1(i,:),'--r',spline_axis,spline_fit.bound2(i,:),'--r','LineWidth',1)
        plot(standard_axis,standard_fit.bound1(i,:),'--g',standard_axis,standard_fit.bound2(i,:),'--g','LineWidth',1)
        title(strcat({'Signal: '},num2str(i)))
        if i == nelectrodes
        h=legend('Spline Model','Standard Model', 'Observed Signal');
        set(h,'Location','SouthEast')
        end
        xlabel('Averaged Spectrum (1/Hz)');
        ylabel('Cumulative Density');
        axis tight
    end
    suptitle('Integrated Spectrum Test')
end

%% Model Test #3 - Test coefficients

if model_true.show_all_plots
    gof_bootstrap_plot( model_true,model_spline,model_standard)
end