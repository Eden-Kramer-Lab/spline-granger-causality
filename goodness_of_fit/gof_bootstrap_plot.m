function gof_bootstrap_plot( model_true,model_spline,model_standard)
% GOF_BOOTSTRAP_PLOT estimates confifence bounds for the coefficient fits using
% a boostrapping procedure. We use the observed coefficient estimates and
% their estimated covariance to generate xxx (nsurrogates variable)
% normally distributedsamples of the coefficients.
% From the resulting distribution, we determine the 0.025 and 0.975
% quantiles of the predicted history dependence.
% In this way, we use the surrogate distribution to define the
% 95% confidence interval for the history dependence estimates.
%
% INPUTS:
% .  model_true     = model structure used to generate the observed data
% .  model_spline   = model fit using spline-Granger method and for which
%                     we compute confidence bounds.
% .  model_standard = model fit using standard-Granger method
%
% OUTPUTS:
% .  plots coefficients for each signal.  Does not plot coefficient fits
%    if neither spline nor standard granger infers an edge and there is no
%    true connection. 

% Initialize
data = model_true.data;
nelectrodes = size(data,1);
model_order = model_true.estimated_model_order;

adj_true= model_true.network;
adj_spline = model_spline.network;
adj_stand = model_standard.network;


b = model_true.true_coefficients;
nlags = size(b,3);  % true model order

bhat = model_spline.model_coefficients;
b_est_stand = model_standard.model_coefficients;

f0 = model_true.sampling_frequency;


dt = 1/f0; % seconds
cntrl_pts = model_true.cntrl_pts;

for electrode = 1:nelectrodes % plot fit for every electrode in network
    
    if sum(adj_spline(electrode,:))==0
        UB = zeros(model_order,nelectrodes);
        LB = zeros(model_order,nelectrodes);
        UB2 = zeros(T*f0,nelectrodes);
        LB2 = zeros(T*f0,nelectrodes);
    else
        [UB,LB]= myBootstrap(model_spline,electrode);
        [ UB2, LB2] = myBootstrapStandard( model_standard, electrode );
    end
    figure;
    k=1;
    for i = 1:nelectrodes
        adj_sum = adj_spline + adj_true + adj_stand;
        adj_sum(adj_sum==2)=1;
        adj_sum(adj_sum==3)=1;
        nconnections = sum(adj_sum(electrode,:));
        if (adj_spline(electrode,i) || adj_true(electrode,i) ||adj_stand(electrode,i))
            subplot(nconnections,1,k)
            plot(dt:dt:(model_order/f0),squeeze(real(bhat(electrode,i,:))),'r','LineWidth',1.5)
            hold on
            plot(dt:dt:(model_order/f0),squeeze(real(b_est_stand(electrode,i,:))),'g','LineWidth',1.5)
            if nlags <= 5
                plot(dt:dt:(nlags/f0),squeeze(real(b(electrode,i,:))),'.k','MarkerSize',30);
            else
                plot(dt:dt:(nlags/f0),squeeze(real(b(electrode,i,:))),'k','LineWidth',2);
            end
            plot(cntrl_pts(2:end)./f0,squeeze(bhat(electrode,i,cntrl_pts(2:end))),'ro','MarkerFaceColor','r')
            plot(dt:dt:(model_order/f0),real(LB(:,i)),'--r','LineWidth',1)
            plot(dt:dt:(model_order/f0),real(UB(:,i)),'--r','LineWidth',1)
          
            plot(dt:dt:(model_order/f0),UB2(:,i),'--g','LineWidth',1.5)
            plot(dt:dt:(model_order/f0),LB2(:,i),'--g','LineWidth',1.5)

            str1 = strcat({'Model Coefficients for Influence of Signal '},num2str(i),{' on Signal '}, num2str(electrode));
            title(str1)
            xlabel('Lag (s)','FontSize',16);
            ylabel('Magnitude','FontSize',16);
            k= k +1;
        end
    end
    h = legend('Spline Estimated','Standard Estimated','True ');
    set(h,'FontSize',15,'Location','NorthEast');
end
end

