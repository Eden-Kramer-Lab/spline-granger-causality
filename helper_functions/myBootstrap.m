function [ UB, LB] = myBootstrap( model, electrode )
% MYBOOTSTRAP creates surrogates for the coefficients estimates
% when building AR models in the spline basis and computes 95% confidence
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
s = model.s;
nsurrogates = model.nsurrogates;
data = model.data;
adj_mat = model.network;
model_order = model.estimated_model_order;
cntrl_pts = model.cntrl_pts;

nelectrodes = size(data,1);            % number electrodes
nobservations = length(data(1,model_order+1:end)); % number of observations

b = zeros(nsurrogates,model_order*nelectrodes);

%%% Define control points and build predictors

c_pt_times = cntrl_pts;  % Define Control Point Locations


% Construct spline regressors for case nelectrodes = 1.
c_pt_times_all = [c_pt_times(1)-100 c_pt_times];
Z = zeros(model_order,length(c_pt_times_all));
num_c_pts = length(c_pt_times_all);  %number of control points in total
for i=1:c_pt_times_all(end-1) %length(t)  %for each 1 ms timepoint, calculate the corresponding row of the glm input matrix
    nearest_c_pt_index = max(find(c_pt_times_all<i));
    nearest_c_pt_time = c_pt_times_all(nearest_c_pt_index);
    next_c_pt_time = c_pt_times_all(nearest_c_pt_index+1);
    
    u = (i-nearest_c_pt_time)/(next_c_pt_time-nearest_c_pt_time);
    
    
    
    p=[u^3 u^2 u 1]*[-s 2-s s-2 s;2*s s-3 3-2*s -s;-s 0 s 0;0 1 0 0];
    Z(i,:) = [zeros(1,nearest_c_pt_index-2) p zeros(1,num_c_pts-4-(nearest_c_pt_index-2))];
    
    
end

for i = c_pt_times_all(end-1)+1:c_pt_times_all(end)
    nearest_c_pt_index = max(find(c_pt_times_all<i));

    nearest_c_pt_time = c_pt_times_all(nearest_c_pt_index);
    next_c_pt_time = c_pt_times_all(nearest_c_pt_index+1);
    
    u = (i-nearest_c_pt_time)/(next_c_pt_time-nearest_c_pt_time);
    
    
    p=[u^3 u^2 u 1]*[-s 2 s-2;2*s -3 3-2*s;-s 0 s;0 1 0];
    Z(i,:) = [zeros(1,nearest_c_pt_index-2) p zeros(1,num_c_pts-4-(nearest_c_pt_index-2))];
    
end



%%% Build model for given electrode
ii = electrode;

% Build history regressors

preds = logical(adj_mat(ii,:)); %% Use results of F-test to input in network
if sum(preds)==0
    bounds = zeros(2,model_order*nelectrodes);
else
    data_copy = data(preds,:); %% Remove electrodes not connected in spline network
    
    
    Z1 = kron(eye(size(data_copy,1)),Z);     % Build block diagonal spline regressors
    
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
    Xfull = X * Z1;      % regressors for y_hat = X*Z1*alpha
    [alpha,~,stats] = glmfit(Xfull,y,'normal','constant','off');  % estimate values at control points, alpha
    
    for kk = 1:nsurrogates
        alpha_hat = alpha + sqrtm(stats.covb)*normrnd(0,1,length(alpha),1);
        % calculate beta values, for every point in space
        bhat(kk,:) = Z1*alpha_hat;                                     % only for electrodes in network
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
    
    max_imaginary_ub = max(imag(bounds(1,:)));
    max_imaginary_lb = max(imag(bounds(2,:)));
    UB = reshape(bounds(1,:),[model_order nelectrodes]);
    LB = reshape(bounds(2,:),[model_order nelectrodes]);
    
end

end

