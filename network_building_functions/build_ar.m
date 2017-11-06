function [ adj_mat] = build_ar( model )
% BUILD_AR builds network model using standard-Granger approach, i.e. 
% condutional Granger-causality
% 
% Outputs:
% . adj_mat = adjacencey matrix for corresponding network


warning off

% 1) Initialize variables and outputs 
data = model.data;
model_order = model.estimated_model_order;
q = model.q;                                       % FDR parameter
nelectrodes = size(data,1);                        % number of electrodes
nobservations = length(data(1,model_order+1:end)); % number of observations
                                                   
                                                  

C = cell(1,nelectrodes);   % cell containing model fit for each electrode
F = zeros(nelectrodes);    % matrix of F-statistics     
X = [];                    % history matrix


% 2) Build history matrix

% Vector of electrode names corresponding to order in build matrix
% e.g. for trivariate case with two lags, [1 1 2 2 3 3]
e_names = zeros(1,model_order*nelectrodes);
j = 1;
for k = 1:nelectrodes
    X_temp = [];
    sgnl = data(k,:)';
    for i=1:model_order                         % For each lag,
        X_temp = [X_temp, circshift(sgnl,i)];   % ... shift x and store it.
        e_names(j) = k;
        j = j+ 1;
    end
    X_temp = X_temp(model_order+1:end,:);
    X = [X X_temp];
end

% 3) Determine connectivity using F test.
%    Build models and test correlation for every electrode pair.

for electrode = 1:nelectrodes
        
    % Generate observations for given y
    x = data(electrode,:);
    y = x(model_order+1:end);
    y = y';
    
    % Fit full model
    [b,~,stats] = glmfit(X,y);
    fit.weights = b;
    fit.pvals = stats.p;
    fit.se = stats.se;
    C{electrode} = fit;
    
    
    % Calculate RSS for full model
    A =[ones(size(X,1),1) X];
    y_hat = A*round(b,10) ;
    error = (y_hat - y).^2;
    rss = sum(error);
    
    % Fit partial model and calculate RSS0 for all minus one electrode
    for ii = 1:nelectrodes
        indices = ~(e_names == ii);
        X0 = X(:,indices);
        [b0,~,stats] = glmfit(X0,y);
        fit0.weights = b;
        fit0.pvals = stats.p;
        fit0.se = stats.se;
        
        A =[ones(size(X0,1),1) X0];
        y_hat = A*round(b0,10);
        error = (y_hat - y).^2;
        rss0 = sum(error);
        
        % Compute F statistic
        F(electrode,ii) = ((rss0 - rss)/model_order)/(rss/(nobservations-nelectrodes*model_order-1));
        
    end
    
end


% 5) Hypothesis test

adj_mat = fpdf(F,model_order,nobservations-nelectrodes*model_order);

m = nelectrodes^2;                 % number of total tests performed
ivals = 1:m;
sp = ivals*q/m;
[pvals, index] = sort(adj_mat(:));
R = find(sp>pvals');               % indices to reject null

adj_mat = zeros(nelectrodes);
adj_mat(index(R)) = 1;             % reject H0 -> correlation

end

