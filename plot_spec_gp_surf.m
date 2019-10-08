% Function to organize and plot data generated from mechanical prediction
% Calls function organize_for_3d_plot to organize data

%function [spec_gp_surf] = plot_spec_gp_surf(metal_ratios,model_number,species_input_model,speciation)

function [spec_gp_surf] = plot_spec_gp_surf(species_input_model, speciation, metal_number)


[~,num_titrations] = size(speciation.titration);

%for metal_number = species_input_model.number.metals + 1
	for titration_number = 1:num_titrations
		% Calls organize_for_3d_plot function for each discrete metal concentration to discretize calculated data for 3D plot
		[spec_gp_mm_3d(titration_number).pH, spec_gp_mm_3d(titration_number).gp] =  organize_for_3d_plot(speciation.titration(titration_number).pH, speciation.titration(titration_number).metal(metal_number).gp_phantom, 2, 13, 0.1);
	end
%end

x = 1:num_titrations;
y = spec_gp_mm_3d(titration_number).pH;

[X,Y] = meshgrid(x,y);

for titration_number = 1:num_titrations
    
    % Makes Z matrix of predicted phantom network moduli
    Z(:,titration_number) = spec_gp_mm_3d(titration_number).gp;
    
end

description_str = char(species_input_model.description);

% Plot 3D Figure
fig_surf = figure;

s = surf(X,Y,Z);
s.EdgeColor = 'interp';
title(description_str,'fontsize',16,'FontName', 'Helvetica Neue')
xlabel('Titration Number','fontsize',14,'FontName', 'Helvetica Neue')
ylabel('pH','fontsize',14,'FontName', 'Helvetica Neue')
zlabel('Plateau Modulus (Pa)','fontsize',14,'FontName', 'Helvetica Neue')

spec_gp_surf.X = X;
spec_gp_surf.Y = Y;
spec_gp_surf.Z = Z;

saveas(fig_surf,description_str)

% Plot Contour Map for nitrocatechol
fontsize = 9.5;


    
    fig_contour_print = figure;
    hold all
    [~,h] = contourf(X,Y,Z./1000,6);
    h.LineWidth = .5;
    h.LineStyle = '-';
    set(gca,'FontName','Helvetica Neue','FontSize',fontsize)
    set(gca,'Box','on')
    xlabel('Titration Number','fontsize',fontsize,'FontName', 'Helvetica Neue','FontWeight','Bold')
    ylabel('pH','fontsize',fontsize,'FontName', 'Helvetica Neue','FontWeight','Bold')
    %xticks([0; 1/3; 2/3; 3/3; 4/3; 5/3]);
    %xticklabels({'0/3'; '1/3'; '2/3'; '3/3'; '4/3'; '5/3'});
    yticks([3; 5; 7; 9; 11; 13]);
    yticklabels({3; 5; 7; 9; 11; 13});
    set(gca,'TickDir','both')
    
    %lim = caxis;
    
%     cmax = 15;
%     cmin = 0;
%     caxis([cmin cmax])
    
    h.LevelList
    c = colorbar;
    c.FontSize = fontsize;
    c.FontName = 'Helvetica Neue';
    c.Label.String = 'Plateau Modulus (kPa)';
    c.Label.FontSize = fontsize;
    c.Label.FontName = 'Helvetica Neue';
    c.Label.FontWeight = 'Bold';
%     c.Ticks = [0  5  10  15];
    c.TickDirection = 'both';
    colormap viridis
    width = 2.75;
    height = 2;
    set(fig_contour_print, 'Position', [400,400,width   *80,height    *76.53])
    
    saveas(fig_contour_print,[description_str, 'contour_print'])
    saveas(fig_contour_print,[description_str, 'contour_print'],'epsc')
    epsclean([description_str, 'contour_print','.eps'])
    
%     if max(h.LevelList) < cmax
%         disp('Well Scaled')
%     else
%         disp('Poorly Scaled')
%     end
    

end