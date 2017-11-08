function plotNetwork( adj)
% PLOTNETWORK plots the adjacency matrix (adj) such that white defines a
%  non-edge and black defines an edge.

imagesc(adj,'AlphaData',0.85); colormap(flipud(gray))

NumTicks = size(adj,2)+1;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
set(gca,'XTickLabel',[])
grid on
names = 1:NumTicks-1;

stringy = {};
for k = 1:length(names)
    stringy = [stringy ; num2str(k)];
end

%set(gca,'XTickLabel',stringy)
L = get(gca,'YLim');
set(gca,'YTick',linspace(L(1),L(2),NumTicks))
set(gca,'YTickLabel',[])
ylabel('node - target', 'FontSize',20);
xlabel('node - source', 'FontSize',20);
axis square

nelectrodes = size(adj,1);
if nelectrodes <= 9
    %%% yticklabel
    xloc = 0.0*ones(1,nelectrodes);
    text(xloc,[1:nelectrodes],stringy,'FontSize',18);
    ylabh = get(gca,'YLabel');
    set(ylabh,'Position',get(ylabh,'Position') - [.5 0 0])
    %%% xticklabel
    yloc = nelectrodes + 0.9*ones(1,nelectrodes);
    text([1:nelectrodes],yloc,stringy,'FontSize',18);
    xlabh = get(gca,'XLabel');
    set(xlabh,'Position',get(xlabh,'Position') + [0 .5 0])
    
end
end

