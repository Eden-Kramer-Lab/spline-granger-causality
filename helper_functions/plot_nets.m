figure;
screensize = get( groot, 'Screensize' );
set(gcf, 'Position', [1,1,screensize(3)/1.5,screensize(4)/3])

subplot 131
plotNetwork(model_true.network);
title('True Configuration','FontSize',16)

subplot 132
plotNetwork(model_standard.network);
title(strcat({'Standard, Runtime: '},num2str(model_standard.computation_time,3),{' s'}),'FontSize',16)

subplot 133
plotNetwork(model_spline.network);
title(strcat({'Spline, Runtime: '},num2str(model_spline.computation_time,3),{' s'}),'FontSize',16)




