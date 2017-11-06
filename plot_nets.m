


if size(model_true.data,1) > 1 % plot if more than one electrode

figure;

subplot 131
plotNetwork(model_true.network);
title('True Configuration','FontSize',16)

subplot 132
plotNetwork(model_standard.network);
title(strcat({'Standard, '},num2str(model_standard.computation_time),{' s'}),'FontSize',16)

subplot 133
plotNetwork(model_spline.network);
title(strcat({'Spline, '},num2str(model_spline.computation_time)),'FontSize',16)

end