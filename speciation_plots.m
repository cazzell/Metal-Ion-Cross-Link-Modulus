% This function is used to plot speciation data generated from calc_speciation and contained in the speciation structure

% This function generates two plots for each of the experimentally relevent metal concentrations
% The first set of plots shows what the polymer bound ligand is coordinating to as a function of pH
% The second set of plots shows what the metal is coordinating to as a function of pH

% This program calls get_spec_plot_structure to organize the data for plotting

function [L_plot,M_plot] = speciation_plots(speciation,species_input_model)

markersize = 6;
markersize_data = 15;
linewidth = 1;
fontsize = 9.5;
titlefontsize = 9.5;
xmin = 1;
xmax = 14;
ymin = 10;
ymax = 10^5;

% Ligand fraction plots
% Ligand row is 1
L_row = 3;

% Get organized data for ligand fraction
[L_plot] = get_spec_plot_structure(L_row,speciation,species_input_model);

for plot_index = 1:length(L_plot.metal_conc)
    
    metal_concentration_str = num2str(speciation.critical_titration.concentrations(plot_index));
    
    fig_L_fraction = figure;
    
    [~,N] = size(L_plot.metal_conc(plot_index).species);
    co = cubehelix(N+2,1.97,-1.45,1.45,1.47,[0.27,0.99],[0.33,1]);
    set(groot,'defaultAxesColorOrder',co);
    
    hold on
    plot(L_plot.metal_conc(plot_index).pH,L_plot.metal_conc(plot_index).species)
    % Uncomment this line to get an overlayed titration curve showing how pH changes as NaOH is added
    % plot(L_plot.metal_conc(plot_index).pH,L_plot.metal_conc(plot_index).norm_v_added)
    lgd = legend(L_plot.metal_conc(plot_index).species_names,'location','best');
    
    xlabel('pH','FontWeight', 'Bold')
    ylabel('Fraction of Ligand','FontWeight', 'Bold')
    
    set(gca,'FontName','Helvetica Neue','FontSize',fontsize)
    set(gca,'Box','on')
    
    legend('boxoff')
    lgd.FontSize = 9.5;
    
    xlimit = xlim;
    ylimit = ylim;
    xmin = xlimit(1);
    xmax = xlimit(2);
    ymin = ylimit(1);
    ymax = ylimit(2);
    axis([0,14,0,1])
    set(gca,'XMinorTick','on')
    set(gca,'YMinorTick','on')
    hA=gca;
    hA.XTick = [0:2:14];
    %hA.XTickLabel = { '0','2','4','6','8','10','12','14'};
    hA.XTickLabel = { '0','2','4','6','8','10','12','14'};
    hA.XAxis.MinorTickValues = [0:1:14];
    hA.YTick = [0:.2:1];
    hA.YTickLabel = { '0.0','0.2','0.4','0.6','0.8','1.0'};
    hA.YAxis.MinorTickValues = [0:.1:1];
    hA.LabelFontSizeMultiplier = 1;
    width = 2.5;
    height = 1.75;
    %set(fig_L_fraction, 'Position', [400,400,width   *80 * 0.9  ,height    * 76.53  *   0.99])
    set(fig_L_fraction, 'Position', [400,400,width   *80 * 0.7  ,height    * 76.53  *   0.99])
	
	axis square
    
    savefig(['speciation_L_frac_',metal_concentration_str,'.fig'])
    saveas(gcf,['speciation_L_frac_',metal_concentration_str,'.eps'],'epsc')
        
end

% Metal fraction plots
% Metal row is 2
M_row = 2;

% Get organized data for metal fraction
[M_plot] = get_spec_plot_structure(M_row,speciation,species_input_model);

for plot_index = 1:length(M_plot.metal_conc)
    
    metal_concentration_str = num2str(speciation.critical_titration.concentrations(plot_index));
    
    fig_M_fraction = figure;
    
    [~,N] = size(M_plot.metal_conc(plot_index).species);
    co = cubehelix(N+2,1.97,-1.45,1.45,1.47,[0.27,0.99],[0.33,1]);
    set(groot,'defaultAxesColorOrder',co);
    
    hold on
    plot(M_plot.metal_conc(plot_index).pH,M_plot.metal_conc(plot_index).species)
    % Uncomment this line to get an overlayed titration curve showing how pH changes as NaOH is added
    % plot(M_plot.metal_conc(plot_index).pH,M_plot.metal_conc(plot_index).norm_v_added)
    lgd = legend(M_plot.metal_conc(plot_index).species_names,'location','best');
    
    xlabel('pH','FontWeight', 'Bold')
    ylabel('Fraction of Metal','FontWeight', 'Bold')
    
    set(gca,'FontName','Helvetica Neue','FontSize',fontsize)
    set(gca,'Box','on')
    
    legend('boxoff')
    lgd.FontSize = 6;
    
    xlimit = xlim;
    ylimit = ylim;
    xmin = xlimit(1);
    xmax = xlimit(2);
    ymin = ylimit(1);
    ymax = ylimit(2);
    axis([0,14,0,1])
    set(gca,'XMinorTick','on')
    set(gca,'YMinorTick','on')
    hA=gca;
    hA.XTick = [0:2:14];
    %hA.XTickLabel = { '0','2','4','6','8','10','12','14'};
    hA.XTickLabel = { '','','4','','8','','12',''};
    hA.XAxis.MinorTickValues = [0:1:14];
    hA.YTick = [0:.2:1];
    hA.YTickLabel = { '0.0','0.2','0.4','0.6','0.8','1.0'};
    hA.YAxis.MinorTickValues = [0:.1:1];
    hA.LabelFontSizeMultiplier = 1;
    width = 2.5;
    height = 1.75;
    %set(fig_M_fraction, 'Position', [400,400,width   *80 * 0.9  ,height    * 76.53  *   0.99])
    set(fig_M_fraction, 'Position', [400,400,width   *80 * 0.7  ,height    * 76.53  *   0.99])
    
    savefig(['speciation_M_frac_',metal_concentration_str,'.fig'])
    saveas(gcf,['speciation_M_frac_',metal_concentration_str,'.eps'],'epsc')
    
end
end
