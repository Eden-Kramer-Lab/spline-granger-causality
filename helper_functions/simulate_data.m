function [model] = simulate_data(model)
% SIMULATE_DATA simulates data using AR model.

b = model.model_coefficients;  % model coefficents
nelectrodes = size(b,1);       % number of electrodes
nlags = size(b,3);             % true model order
T = model.T;                   % total length of recording (seconds)
f0 = model.sampling_frequency;


N = T*f0 + nlags + 3000;      % simulate 3000 additional points
noise = model.noise;
X = noise.*randn([nelectrodes N]);

for t = nlags+1:N
    for k = 1:nlags
        X(:,t) = X(:,t) + b(:,:,k)*X(:,t-k);
    end
end
X = X(:,nlags+3001:end);      % remove first 3000
model.data = X;

end
