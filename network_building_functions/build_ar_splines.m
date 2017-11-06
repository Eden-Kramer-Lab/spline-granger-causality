function [ adj_mat] = build_ar_splines( model )
% BUILD_AR_SPLINES builds network model using spline-Granger approach and
% outputs:
% . adj_mat = adjacencey matrix for corresponding network
%

warning off

% 1) Initialize variables & outputs 
data = model.data;
model_order = model.estimated_model_order;
cntrl_pts = model.cntrl_pts;
q = model.q; % max number acceptable proportion of false discoveries
s = model.s; % tension parameter
nobservations = length(data(1,model_order+1:end)); % number of observations
nelectrodes = size(data,1);                        % number of electrodes
F = zeros(nelectrodes);                            % matrix of F-statistics

% 2) Define control points and build predictors
%    Construct spline regressors.

c_pt_times_all = [cntrl_pts(1)-100 cntrl_pts];
Z = zeros(model_order,length(c_pt_times_all));
num_c_pts = length(c_pt_times_all);     % number of control points in total
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
Z0 = kron(eye(nelectrodes-1),Z);   % Nested model spline regressors
Z1 = kron(eye(nelectrodes),Z);     % Full model spline regressors

% 3) Build history matrix
X = [];

% Vector of electrode names corresponding to order in build matrix
% e.g. for trivariate case with two lags, [1 1 2 2 3 3]
e_names = zeros(1,model_order*nelectrodes);
j = 1;

for k = 1:nelectrodes
    X_temp = [];
    sgnl = data(k,:)';
    for i=1:model_order                         % For each lag,
        X_temp = [X_temp, circshift(sgnl,i)];   %... shift x and store it.
        e_names(j) = k;
        j = j+ 1;
    end
    X_temp = X_temp(model_order+1:end,:);
    X = [X X_temp];
end

% 4) Determine connectivity using F test.
%    Build models and test correlation for every electrode pair.

for electrode = 1:nelectrodes
    
    % Generate observations for given y
    x = data(electrode,:);
    y = x(model_order+1:end);
    y = y';
    
    % Fit full model and calculate RSS
    Xfull = X * Z1;      % regressors for y_hat = X*Z1*alpha
    [alpha,~,stats] = glmfit(Xfull,y,'normal','constant','off');
    fit.weights = alpha;
    fit.pvals = stats.p;
    fit.se = stats.se;
    A = Xfull;
    y_hat = A*round(alpha,10);
    error = (y_hat - y).^2;
    rss = sum(error);
    
    % Fit partial model and calculate RSS0 for all minus one electrode
    for ii = 1:nelectrodes
        indices = ~(e_names == ii);
        X0 = X(:,indices); 
        X0full = X0 * Z0;  
        [a0,~,stats] = glmfit(X0full,y,'normal','constant','off');
        fit0.weights = a0;
        fit0.pvals = stats.p;
        fit0.se = stats.se;
        A = X0full;
        y_hat = A*round(a0,10) ;
        error = (y_hat - y).^2;
        rss0 = sum(error);
        
        % Compute F statistic
        F(electrode,ii) = ((rss0 - rss)/num_c_pts)/(rss/(nobservations-nelectrodes*num_c_pts));
        
    end
    
end


% 5) Hypothesis test

adj_mat = fpdf(F,num_c_pts,nobservations-nelectrodes*num_c_pts);
m = nelectrodes^2;                  % total number of tests performed
ivals = 1:m;
sp = ivals*q/m;
[pvals, index] = sort(adj_mat(:));
R = find(sp>pvals');               % indices to reject null
adj_mat = zeros(nelectrodes);
adj_mat(index(R)) = 1;             % reject H0 -> correlation

end

