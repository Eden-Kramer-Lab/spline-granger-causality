function [notwhite, dw] = gof_residuals( model)
% GOF_RESIDUALS analyzes the goodness of fit of the model residuals. For
% each signal, it plots the true signal and the estimated signal, the model
% residuals over time and an autocorrelation plot of the residuals.
%
% INPUTS:
% . model = structure containing data and estimated data.
%
% OUTPUTS:
% . notwhite = signals that fail the Durbin-Watson test for serial
%              autocorrelation, implying poor model fit.
% . dw       = the Durbin-Watson statistic

% Define and initialize parameters
data  = model.data;
nelectrodes = size(data,1);
model_order = model.estimated_model_order;
datap = data(:,model_order+1:end);
yestimate = model.signal_estimate;
dt = 1/model.sampling_frequency;
T = model.T;
taxis = dt:dt:T';
residuals = zeros(nelectrodes,length(datap));

% Analyze residuals for each electrode fit
for electrode = 1:nelectrodes
    
    true_signal =  datap(electrode,:);   
    est_signal = yestimate(electrode,:);
    residuals(electrode,:) = true_signal-est_signal;
    
    if model.show_all_plots
        h= figure;
        % Plot the true and estimated/recondstructed signal
        subplot 311
        plot(taxis(model_order+1:end),true_signal,'k');
        hold on;
        plot(taxis(model_order+1:end),est_signal,'r');
        xlabel('Time(s)');
        legend('true signal','estimated signal')
        
        % Calculuate and plot residuals
        
        subplot 312
        plot(taxis(model_order+1:end),residuals(electrode,:),'.');
        plot(residuals(electrode,:),'.');
        title('Residuals')
        xlabel('Time(s)');
        
        % Plot autocorrelation of residuals
        subplot 313
        autocorr(residuals(electrode,:));
        title(['Fit for Signal: ' num2str(electrode) ' of ' model.name])
    end
end

% Calculate DW- statistic
[dw, pval] = whiteness(data,residuals);

% Control for mulitple hypothesis using False-Discovery Rate
sig = significance(pval,0.05,'FDR');
notwhite = find(sig);

end

