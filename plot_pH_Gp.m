
% for metal_concentration = 1/3;
% I ran this with model_2 speciation data loaded to generate computational
% schematic figure

index = speciation.critical_titration.indices(1);
pH = speciation.pH;
gp = speciation.metal_concentration(index).gp_phantom;


fig_gp_pH = figure;

hold on
plot(pH,gp./1000)

xlabel('pH','FontWeight', 'Bold')
ylabel('Modulus','FontWeight', 'Bold')

set(gca,'FontName','Helvetica Neue','FontSize',7)
set(gca,'Box','on')

legend('boxoff')
lgd.FontSize = 7;

xlimit = xlim;
ylimit = ylim;
xmin = xlimit(1);
xmax = xlimit(2);
ymin = ylimit(1);
ymax = ylimit(2);
axis([0,14,0,18])
set(gca,'XMinorTick','on')
set(gca,'YMinorTick','on')
hA=gca;
hA.XTick = [0:2:14];
%hA.XTickLabel = { '0','2','4','6','8','10','12','14'};
hA.XTickLabel = { '0','2','4','6','8','10','12','14'};
hA.XAxis.MinorTickValues = [0:1:14];
hA.YTick = [0:5:18];
hA.YTickLabel = { '0','5','10','15'};
hA.YAxis.MinorTickValues = [0:1:18];
hA.LabelFontSizeMultiplier = 1;

width = 2.5;
height = 1.75;
set(fig_gp_pH, 'Position', [400,400,width   *80 * 0.7  ,height    * 76.53  *   0.99])

axis square

savefig(['gp_pH','.fig'])
saveas(gcf,['gp_pH','.eps'],'epsc')

