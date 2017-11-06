
%%% Simulation parameters -------------------------------------------------

model_true.sampling_frequency = 500;
model_true.T = 2;   % time in seconds of window
model_true.noise = 0.25; % standard deviation of the noise
model_true.taxis = (1/model_true.sampling_frequency):(1/model_true.sampling_frequency):model_true.T;
model_true.true_coefficients = single_node_order20;
model_true.model_coefficients = model_true.true_coefficients;   

%%% Define model inputs for spline Granger & standard Granger -------------

model_true.s = 0.5;                     % spline tenstion paramter
model_true.estimated_model_order = 20;  % history dependence in model 
                                        % ... (samples)
model_true.cntrl_pts =0:5:model_true.estimated_model_order; % control points

%%% Define network testing parameters -------------------------------------

model_true.q = 0.05;            % FDR parameter: acceptable proportion of 
                                % ... false discoveries


