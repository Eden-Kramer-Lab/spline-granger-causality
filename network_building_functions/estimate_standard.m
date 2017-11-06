function [bhat, yhat] = estimate_standard(model, adj_mat )
% ESTIMATE_STANDARD fits standard-Granger model given an predetermined
% network structure.
%
% INPUTS:
%  model           = model structure containing data (with dimensions
%                    electrodes x time), and estimated model order.
%  adj_mat        =  desired network configuration

%
% OUTPUTS:
%  bhat   = coefficient estimates
%  yhat    = signal estimates

% Initialize variables & outputs
data = model.data;
model_order = model.estimated_model_order;
nelectrodes = size(data,1);                        % number electrodes
nobservations = length(data(1,model_order+1:end)); % number of observations
bhat = zeros(nelectrodes,nelectrodes,model_order);
yhat = zeros(nelectrodes,size(data,2));
yhat = yhat(:,model_order+1:end);

% Fit standard-Granger MVAR model for each electrode
for electrode = 1:nelectrodes
    model_size = sum(adj_mat(electrode,:));
    if model_size > 0
        % Generate observations for given y
        x = data(electrode,:);
        y = x(model_order+1:end);
        y = y';

        % Build regressors for specified network fit.  
        indices = logical(adj_mat(electrode,:)); % Find connected nodes,
        data_subset = data(indices,:);           % ... and only use data for
                                                 % ... those nodes
        X = []; % History matrix
        for k = 1:model_size
            X_temp = [];
            sgnl = data_subset(k,:)';
            for i=1:model_order                         % For each lag,
                X_temp = [X_temp, circshift(sgnl,i)];   %... shift x and store it.
            end
            X_temp = X_temp(model_order+1:end,:);
            X = [X X_temp];
        end
        [b,~,~] = glmfit(X,y,'normal','constant','off'); % Fit model
        yhat(electrode,:) = glmval(b,X,'identity','constant','off'); % Get signal estimate.
        
        % Format inferred coefficients appropriately
        j=1;
        for p = 1:nelectrodes
            if adj_mat(electrode,p) == 1
                bhat(electrode,p,:) = b(j:j+model_order-1);
                j= j+model_order;
            end
        end
    end 
end
end

