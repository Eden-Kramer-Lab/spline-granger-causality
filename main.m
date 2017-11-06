clear all;
close all;

%%% Configure simulation --------------------------------------------------
config;

%%% Simulate data ---------------------------------------------------------
model_true = simulate_data(model_true);

%%% Infer networks using standard-Granger and spline-Granger---------------

[model_true, model_spline, model_standard] = infer_nets(model_true);

%%% Plot results ----------------------------------------------------------
plot_nets;

%%% Check goodness of fit -------------------------------------------------
model_testing;

%%%%%% QUESTIONS
% x. Include bootstrap coefficients? yes
% x. Include plots of autocorrelation plots? yes
% x. whiteness, significance and DEmean are functions from MVGC toolbox add
% note
