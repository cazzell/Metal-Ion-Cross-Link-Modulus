% Function to organize and plot data generated from mechanical prediction
% Calls function organize_for_3d_plot to organize data

function [spec_gp_surf] = plot_spec_gp_surf_carrr(metal_ratios,species_input_model,speciation,mechanics,metal_number,ligand_number)

[~, num_titrations] = size(speciation.titration_number);

x = 1:num_titrations;
y = species_input_model.pH;

[X,Y] = meshgrid(x,y);

for titration_number = 1:length(num_titrations)
	% Makes Z matrix of predicted phantom network moduli
	Z(:,titration_number) = mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).gp_phantom;
end

description_str = char(species_input_model.description);

spec_gp_surf.X = X;
spec_gp_surf.Y = Y;
spec_gp_surf.Z = Z;

% Plot 3D Figure
fig_surf = figure;

s = surf(X,Y,Z);
s.EdgeColor = 'interp';
title(description_str,'fontsize',16,'FontName', 'Helvetica Neue')
xlabel('Metal Equivalency','fontsize',14,'FontName', 'Helvetica Neue')
ylabel('pH','fontsize',14,'FontName', 'Helvetica Neue')
zlabel('Plateau Modulus (Pa)','fontsize',14,'FontName', 'Helvetica Neue')

axis([0,5/3,0,14,0,17000])
xticks([0; 1/3; 2/3; 3/3; 4/3; 5/3]);
xticklabels({'0/3'; '1/3'; '2/3'; '3/3'; '4/3'; '5/3'});
yticks([0; 2; 4; 6; 8; 10; 12; 14]);
yticklabels({0; 2; 4; 6; 8; 10; 12; 14});
hA = gca;
hA.YAxis.MinorTickValues = [0:1:14];

axis square

set(gcf,'renderer','Painters')

saveas(fig_surf,description_str)

    saveas(fig_surf,[description_str, 'surf_print'])
    saveas(fig_surf,[description_str, 'surf_print'],'epsc')

	fontsize = 9.5;
    
    fig_contour_print = figure;
    hold all
    [~,h] = contourf(X,Y,Z./1000,6);
    h.LineWidth = .5;
    h.LineStyle = '-';
    set(gca,'FontName','Helvetica Neue','FontSize',fontsize)
    set(gca,'Box','on')
    xlabel('Metal Equivalency','fontsize',fontsize,'FontName', 'Helvetica Neue','FontWeight','Bold')
    ylabel('pH','fontsize',fontsize,'FontName', 'Helvetica Neue','FontWeight','Bold')
    xticks([0; 1/3; 2/3; 3/3; 4/3; 5/3]);
    xticklabels({'0/3'; '1/3'; '2/3'; '3/3'; '4/3'; '5/3'});
    yticks([0; 2; 4; 6; 8; 10; 12; 14]);
    yticklabels({0; 2; 4; 6; 8; 10; 12; 14});
    %set(gca,'TickDir','both')
	hA = gca;
	hA.YAxis.MinorTickValues = [0:1:14];
    
    %lim = caxis;
    
    cmax = 15;
    cmin = 0;
    caxis([cmin cmax])
    
    h.LevelList
    c = colorbar;
    c.FontSize = fontsize;
    c.FontName = 'Helvetica Neue';
    c.Label.String = 'Plateau Modulus (kPa)';
    c.Label.FontSize = fontsize;
    c.Label.FontName = 'Helvetica Neue';
    c.Label.FontWeight = 'Bold';
    c.Ticks = [0  5  10  15];
    c.TickDirection = 'both';
    %colormap viridis
    width = 2.75;
    height = 2;
    set(fig_contour_print, 'Position', [400,400,width   *80,height    *76.53])
    
	axis square
	
	set(gcf,'renderer','Painters')
	
    saveas(fig_contour_print,[description_str, 'contour_print'])
    saveas(fig_contour_print,[description_str, 'contour_print'],'epsc')
    
    if max(h.LevelList) < cmax
        disp('Well Scaled')
    else
        disp('Poorly Scaled')
    end
    
end