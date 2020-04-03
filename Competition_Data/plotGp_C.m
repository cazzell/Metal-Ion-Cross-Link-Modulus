
% M_1_C_03_mech = load('M_1_C_03_mech.mat');
% M_2_C_03_mech = load('M_2_C_03_mech.mat');
% M_3_C_03_mech = load('M_3_C_03_mech.mat');

M_1_C_03_mech = load('M_1_C_02_mech.mat');
M_2_C_03_mech = load('M_2_C_02_mech.mat');
M_3_C_03_mech = load('M_3_C_02_mech.mat');

% pH = 7 at 126
pH_number = 126;

x = linspace(0,2,100);
M_1_y = M_1_C_03_mech.mechanics.plot_gp.ligand(3).metal(2).Z(pH_number,:);
M_2_y = M_2_C_03_mech.mechanics.plot_gp.ligand(3).metal(2).Z(pH_number,:);
M_3_y = M_3_C_03_mech.mechanics.plot_gp.ligand(3).metal(2).Z(pH_number,:);

figure
hold on
plot(x,M_1_y./1000)
plot(x,M_2_y./1000)
plot(x,M_3_y./1000)

axis([0,2,0,25])
set(gca,'Box','on')

axis square