function [ LB, UB] = myBootstrapStandard( model, electrode )
% MYBOOTSTRAP creates surrogates for the coefficients estimates
% when building AR models in the standard basis and computes 95% confidence
% intervals.
%
% INPUTS:
%  1. electrode  = the electrode whose data is used for the model fit
%  2. model      = structure that contains the following:
%      data         = A matrix of electode data with dimensions electrodes
%                     x time
%      adj_mat      = adjacencey matrix for corresponding network
%      model_order  = The number of lags used as used for predictor variables
%      nsurrogates  = number of surrogates
%
% OUTPUT
%  UB: upper bounds
%  LB: lower bounds
nsurrogates = model.nsurrogates;
data = model.data;
adj_mat = model.network;
model_order = model.estimated_model_order;


nelectrodes = size(data,1);            % number electrodes

b = zeros(nsurrogates,model_order*nelectrodes);

%%% Define control points and build predictors

%%% Build model for given electrode
ii = electrode;

% Build history regressors
preds = logical(adj_mat(ii,:)); %% Use results of F-test to input in network
if sum(preds)==0
    bounds = zeros(2,model_order*nelectrodes);
     LB = reshape(bounds(1,:),[model_order nelectrodes]);
   UB = reshape(bounds(2,:),[model_order nelectrodes]);
else
    data_copy = data(preds,:); %% Remove electrodes not connected in spline network
    X = [];                                 % Build history matrix
    for k = 1:size(data_copy,1)
        X_temp = [];
        sgnl = data_copy(k,:)';
        for i=1:model_order                                   %For each lag,
            X_temp = [X_temp, circshift(sgnl,i)];   %... shift x and store it.
        end
        X_temp = X_temp(model_order+1:end,:);
        X = [X X_temp];
    end
    
    
    %%% Build Model
    
    % Generate observations for given y
    y = data(ii,model_order+1:end);
    y = y';
    
    % Fit full model and calculate RSS
    Xfull = X;      % regressors for y_hat = X*Z1*alpha
    %[alpha,~,stats] = glmfit(Xfull,y,'normal','constant','off');  % estimate values at control points, alpha
        [mdl1] = fitglm(Xfull,y,'Distribution','normal','Intercept',false);
        
        alpha = mdl1.Coefficients.Estimate;
        covb = mdl1.CoefficientCovariance;
    for kk = 1:nsurrogates
        alpha_hat = alpha + sqrtm(covb)*normrnd(0,1,length(alpha),1);
        % calculate beta values, for every point in space
        bhat(kk,:) = alpha_hat;                                     % only for electrodes in network
        j =1;
        for k = 1:nelectrodes
            if preds(k)
                b(kk,((k-1)*model_order + 1): k*model_order) = bhat(kk,((j-1)*model_order + 1): j*model_order);
                j = j+1;
            end
        end
        
    end
    
    ind1 =round(nsurrogates*0.025);
    ind2 = round(nsurrogates*0.975);
    
    if ind1==0
        ind1=1;
    end
    sorted_b = sort(real(b));
    bounds(1,:) = sorted_b(ind1,:);
    bounds(2,:) = sorted_b(ind2,:);
    LB = reshape(bounds(1,:),[model_order nelectrodes]);
    UB = reshape(bounds(2,:),[model_order nelectrodes]);
    
end

end

