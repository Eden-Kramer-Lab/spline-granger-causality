%%%%------------------ Configure Simulation--------------------------------

model_true.show_all_plots = true;       % if true, all plots for all signal fits will be
                                        % ... will be output.  If false, only print-out of
                                        % ... results will be shown.
                       
%%% Simulation parameters
model_true.true_coefficients = single_node_high_freq ;
model_true.model_coefficients = model_true.true_coefficients;   
model_true.sampling_frequency = 500;    % sampling frequency of signal in Hertz.
model_true.T = 2;                       % time in seconds of window
model_true.noise = 0.25;                % standard deviation of the noise
model_true.taxis = (1/model_true.sampling_frequency):(1/model_true.sampling_frequency):model_true.T;

%%% Define model inputs for spline Granger & standard Granger

model_true.s = 0.5;                     % spline tension paramter
model_true.estimated_model_order = 30;  % history dependence in model (samples)
model_true.cntrl_pts =0:5:model_true.estimated_model_order; % control points

%%% Define network testing parameters

model_true.q = 0.05;                    % FDR parameter: acceptable proportion of false discoveries
model_true.nsurrogates = 10000;         % Number of surrogates used for bootstrap coefficients
                    
%%%---------------------- Simulate data -----------------------------------
model_true = simulate_data(model_true);

%%%-------- Infer networks using standard-Granger and spline-Granger-------

[model_true, model_spline, model_standard] = infer_nets(model_true);

%%% ---------------------- Check goodness of fit -------------------------
model_testing;


