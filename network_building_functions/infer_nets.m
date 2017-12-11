
function [model_true, model_spline, model_standard] = infer_nets(model_true)
%%% Determine true network connectivity -----------------------------------

adj_true = model_true.true_coefficients;
adj_true(adj_true~=0)=1;
adj_true=sum(adj_true,3);
adj_true(adj_true~=0)=1;
model_true.network = adj_true;

%%% Fit spline-Granger model to data --------------------------------------
tic
[ adj_spline] = build_ar_splines( model_true);
splinetime  = toc;
[ bhat, yhat] = estimate_coefficient_fits( model_true, adj_spline);

model_spline = model_true;
model_spline.model_coefficients = bhat;
model_spline.computation_time = splinetime;
model_spline.signal_estimate = yhat;
model_spline.network = adj_spline;
model_spline.name = 'Spline-Granger Model';


%%% Fit standard-Granger model to data ------------------------------------
tic
[ adj_standard] = build_ar( model_true );
standardtime  = toc;
[ bhat, yhat] = estimate_standard( model_true, adj_standard);

model_standard = model_true;
model_standard.model_coefficients = bhat;
model_standard.computation_time = standardtime;
model_standard.signal_estimate = yhat;
model_standard.network = adj_standard;
model_standard.name = 'Standard-Granger Model';


end


