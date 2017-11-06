function [ model_fit] = grstat( model_true,model_estimated)
% GRSTAT implements the Grenander and Rosenblatt test of the integrated
% spectrum, [REF = Spectral Analysis and Time Series, Volume 1: Univariate
% Series by M. B. Priestley].
%
% INPUTS:
% . model_true      = the model from which the observed data was generated
% . model_estimated = the model, either spline-Granger or standard-Granger,
%                     that was fit to the observed data
%
% OUTPUTS:
% . model_fit = structure containing assessment of the fit model,
%               model_estimated's fit to the observed data, model_true.
%               Structure containing the integrated spectrum of both models
%               and 95% confidence bounds for the fit signal. Contains
%               the GR-statistic, and associated pvalue.  Contains signals
%               that fail the test (implying poor model fit) using
%               Bonferroni correction.


% Define and intialize parameters
nelectrodes = size(model_true.data,1); % number of electrodes
f0 = model_true.sampling_frequency;
trials_true = model_true.data;

% Generate new realization of process from estimated model
model_estimated_new = simulate_data(model_estimated);
trials_estimated = model_estimated_new.data;

% Analyze spectra for each electrode fit
for electrode = 1:nelectrodes
    x=trials_true(electrode,:);
    x = x - ones(size(x)).*mean(x);
    TW = 2; % time-bandwidth product
    [Sxx,~] = pmtm(x,TW,length(x),f0); % Compute spectrum for observed signal
    h_true = Sxx;
    
    
    x=trials_estimated(electrode,:); % Compute spectrum for estimated model signal
    x = x - ones(size(x)).*mean(x);
    [Sxx, ~] = pmtm(x,TW,length(x),f0);
    h_estimated  = Sxx;
    
    [H, X] = ecdf(h_true);           % Compute CDF of spectra for true
    [H1, X1] = ecdf(h_estimated);    % ... and for estimated model signal

    ap = 2.2414;                     % Defines 95% confidence
    total_observations = size(model_true.data,2);  
    
    flag = 'biased';                              % divide by 1/N
    R = xcov(trials_estimated(electrode,:),flag); % autocovariance of estimated signal
    
    G = sum(R(970:end-30).^2);
    G = G/(4*pi);
    
    conf1 = ap*sqrt(8*pi*G/total_observations);   % Define confidence bounds
    UB = H1 + conf1;                              % ... Upper bound
    LB = H1 - conf1;                              % ... Lower bound
    
    % Linear interpolation to fill in missing values
    q1 = [X' X1'];
    Q = NaN(5,length(q1));
    Q(1,:) = q1;
    Q(2,:) = interp1q(X,H,Q(1,:)')';
    Q(3,:) = interp1q(X1,LB,Q(1,:)')';
    Q(4,:) = interp1q(X1,H1,Q(1,:)')';
    Q(5,:) = interp1q(X1,UB,Q(1,:)')'';
    
    Q = sortrows(Q',1)';
    Qp=Q';
    Qp = Qp(~any(isnan(Qp),2),:)';
    
    % Compute GR-statistic
    gr_estimated = max(sqrt(total_observations)*abs(Qp(2,:)-Qp(4,:)));
    
    % Store values
    model_fit.xaxis(electrode,:) = Q(1,:);
    model_fit.true(electrode,:) = Q(2,:);
    model_fit.estimate(electrode,:) = Q(4,:);
    model_fit.bound1(electrode,:) = Q(3,:);
    model_fit.bound2(electrode,:) = Q(5,:);
    
    if isnan(gr_estimated)
        gr_estimated=0;
    end
    model_fit.stat(electrode) = gr_estimated;
    
end

% Control for problems of multiple testing using Bonferonni correction
bfcorrection = 0.05/length(model_fit.stat);
model_fit.pvals =  pval2grstat(model_fit.stat,'grstatistic');
model_fit.fails = find(model_fit.pvals<bfcorrection);
end

