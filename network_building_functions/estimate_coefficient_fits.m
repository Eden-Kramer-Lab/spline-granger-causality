function [bhat, yhat] = estimate_coefficient_fits( model, adj_mat)
% ESTIMATE_COEFFICIENT_FITS fits spline-Granger model given an predetermined
% network structure.
%
% INPUTS:
%  model           = model structure containing data (with dimensions
%                    electrodes x time), estimated model order, and the
%                    control points.
%  adj_mat        =  desired network configuration

%
% OUTPUTS:
%  bhat   = coefficient estimates
%  yhat    = signal estimates


% Initialize variables & outputs
s = model.s;                                       % tension parameter
data= model.data;
model_order = model.estimated_model_order;
cntrl_pts  = model.cntrl_pts;
nelectrodes = size(data,1);                        % number electrodes
nobservations = length(data(1,model_order+1:end)); % number of observations
bhat = zeros(nelectrodes,nelectrodes,model_order);
yhat = zeros(nelectrodes,size(data,2));
yhat = yhat(:,model_order+1:end);

% Define transformation matrix into spline basis
c_pt_times_all = [cntrl_pts(1)-100 cntrl_pts];
Z = zeros(model_order,length(c_pt_times_all));
num_c_pts = length(c_pt_times_all);  % number of control points in total
for i=1:c_pt_times_all(end-1)        
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

% Fit spline-Granger MVAR model for each electrode
for electrode = 1:nelectrodes
    model_size = sum(adj_mat(electrode,:));
    if model_size > 0
        Z0 = kron(eye(model_size),Z);   % Full matrix with spline regressors
        
        % Generate observations for given y
        x = data(electrode,:);
        y = x(model_order+1:end);
        y = y';

        % Build regressors for specified network fit.
        indices = logical(adj_mat(electrode,:));
        data_subset = data(indices,:);
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
        Xfull = X * Z0;  % Regressors in spline basis
        [alpha,~,stats] = glmfit(Xfull,y,'normal','constant','off'); % Fit model
        yhat(electrode,:) = glmval(alpha,Xfull,'identity','constant','off'); % Get signal estimate.
        b = Z0*alpha;    % Transform out of spline basis
        
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

