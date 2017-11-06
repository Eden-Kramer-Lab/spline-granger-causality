
function b = single_node_low_freq
% Returns coefficients for univariate AR(1) or AR(2) model that generate a
% low frequency signal
%

model_coefficients = [0.3 0.3];    % AR(2) coefficients
% model_coefficients = [0.9 -0.1]; % ...alternate AR(2) coefficients
% model_coefficients = 0.9;        % ...alternate AR(1) coefficients


b = zeros(1,1,length(model_coefficients));
b(1,1,:) = model_coefficients;

end