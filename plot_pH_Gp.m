%% Code authored by Seth Allen Cazzell
% cazzell.lbi@gmail.com

% This is the program to plot the mechanical prediction at a fixed pH.

%% Please reference
% "Expanding the stoichiometric window for metal cross-linked gel assembly using competition" PNAS, Cazzell (2019)
% when appropriate.

%% Copyright 2019 Seth Allen Cazzell All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%%

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

