%% Code authored by Seth Allen Cazzell
% cazzell.lbi@gmail.com

% This is the program to plot the mechanical predictions. This function is
% called by "calc_mechanics"

%% Please reference
% "Expanding the stoichiometric window for metal cross-linked gel assembly using competition" PNAS, Cazzell (2019)
% when appropriate.

%% Copyright 2019 Seth Allen Cazzell All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%%
function [X,Y,Z] = plot_gp_contour(species_input_model,speciation,mechanics,metal_number,ligand_number)

[~, num_titrations] = size(speciation.titration_number);

x = 1:num_titrations;
y = speciation.pH_list;

[X,Y] = meshgrid(x,y);

for titration_number = 1:num_titrations
	% Makes Z matrix of predicted phantom network moduli
	Z(:,titration_number) = mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).gp_phantom;
end

description_str = char(species_input_model.description);

% plot_gp.ligand(ligand_number).metal(metal_number).X = X;
% plot_gp.ligand(ligand_number).metal(metal_number).Y = Y;
% plot_gp.ligand(ligand_number).metal(metal_number).Z = Z;

% % Plot 3D Figure
% fig_surf = figure;
%
% s = surf(X,Y,Z);
% s.EdgeColor = 'interp';
% title(description_str,'fontsize',16,'FontName', 'Helvetica Neue')
% xlabel('Metal Equivalency','fontsize',14,'FontName', 'Helvetica Neue')
% ylabel('pH','fontsize',14,'FontName', 'Helvetica Neue')
% zlabel('Plateau Modulus (Pa)','fontsize',14,'FontName', 'Helvetica Neue')
%
% axis([0,5/3,0,14,0,17000])
% xticks([0; 1/3; 2/3; 3/3; 4/3; 5/3]);
% xticklabels({'0/3'; '1/3'; '2/3'; '3/3'; '4/3'; '5/3'});
% yticks([0; 2; 4; 6; 8; 10; 12; 14]);
% yticklabels({0; 2; 4; 6; 8; 10; 12; 14});
% hA = gca;
% hA.YAxis.MinorTickValues = [0:1:14];
%
% axis square
%
% set(gcf,'renderer','Painters')
%
% saveas(fig_surf,description_str)
%
% saveas(fig_surf,[description_str, 'surf_print'])
% saveas(fig_surf,[description_str, 'surf_print'],'epsc')

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
yticks([0; 2; 4; 6; 8; 10; 12; 14]);
yticklabels({0; 2; 4; 6; 8; 10; 12; 14});
%set(gca,'TickDir','both')
hA = gca;
hA.YAxis.MinorTickValues = [0:1:14];

%lim = caxis;

% cmax = 15;
% cmin = 0;
% caxis([cmin cmax])

h.LevelList;
c = colorbar;
c.FontSize = fontsize;
c.FontName = 'Helvetica Neue';
c.Label.String = 'Plateau Modulus (kPa)';
c.Label.FontSize = fontsize;
c.Label.FontName = 'Helvetica Neue';
c.Label.FontWeight = 'Bold';
%c.Ticks = [0  5  10  15];
c.TickDirection = 'both';
%colormap viridis
width = 2.75;
height = 2;
set(fig_contour_print, 'Position', [400,400,width   *80,height    *76.53])

axis square
set(gcf,'renderer','Painters')

saveas(fig_contour_print,[description_str, 'contour_print'])
saveas(fig_contour_print,[description_str, 'contour_print'],'epsc')

% if max(h.LevelList) < cmax
% 	disp('Well Scaled')
% else
% 	disp('Poorly Scaled')
% end

end