%% Model Test #1 - Durbin-Watson test for autocorrelation of model residuals

gof_residuals(model_spline);
title('Spline Fit')
gof_residuals(model_standard);
title('Standard Fit')


%% Model Test #2 - Grenander and Rosenblatt test of the Integrated Spectrum

[spline_fit] = grstat(model_true,model_spline);
[standard_fit] = grstat(model_true,model_standard);
