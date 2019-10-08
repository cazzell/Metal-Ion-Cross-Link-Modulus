% Code authored by Seth Cazzell
% cazzell.lbi@gmail.com

% Code used to calculate the speciation of metal-coordinate complexes as a function of both pH and metal concentration
% to predict the plateau moduli as a function of these variables

% This is the master program. Here the variable to modify is the model_number.

close all
clear
clc

warning('off','all')

% Generates model used to describe the desired system for the speciation program
[species_input_model] = generate_models()

% Runs speciation program to predict species concentrations vs. pH and metal concentration
[speciation] = calc_speciation(species_input_model)

% % Function to organize and generate plots
% for metal_number = 1:species_input_model.number.metals+1
% 	[spec_gp_surf] = plot_spec_gp_surf(species_input_model, speciation, metal_number)
% end


% Generates plots to show both the fraction of what the ligand is bound to and the fraction of what the metal is bound to
%[L_spec_plot,M_spec_plot] = speciation_plots(speciation,species_input_model)

% Function to organize and generate plots
%[spec_gp_surf] = plot_spec_gp_surf(metal_ratios,model_number,species_input_model,speciation)

% Code to save generated data
save('structures.mat','species_input_model','speciation')
%save('structures.mat','species_input_model','speciation','spec_gp_surf')


