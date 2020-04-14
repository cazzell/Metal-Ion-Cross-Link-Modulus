%% Code authored by Seth Allen Cazzell
% cazzell.lbi@gmail.com

% This is the program to plot speciation data from already calculated data
% First run model to perform calculation in gel_mechanics_model.m

%% Please reference
% "Expanding the stoichiometric window for metal cross-linked gel assembly using competition" PNAS, Cazzell (2019)
% when appropriate.

%% Copyright 2020 Seth Allen Cazzell All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% Setup and Clear Space
close all
clear
clc

%% Things you can plot
% First, this currently just works if you want to plot systems built with one ligand and one metal

% From the contour plots there are two variables, y-axis is pH
% x-axis is what you titrate, like metal or ligand concentration,
% or some other compositionial change. I'll refer to this as the titration variable

% So we can plot modulus or speciation with either a constant pH or
% constant titration variable. To give a total of 4 different styles of plots

% Modulus as a function of pH and constant titration variable
% Modulus as a function of titration variable (x-axis of contour plot) and constant pH

% Speciation as a function of pH and constant titration variable
% Speciation as a function of titration variable (x-axis of contour plot) and constant pH

%% Inputs
% Define y-axis for this plot to be either modulus or speciation
% speciation is 1
% modulus is 2
y_choice = 2;

% Constant pH, or other titration variable
% constant pH is 1
% constant titration is 2
constant_choice = 2;

% regardless of choice, specify the component number that you want to examine
component_number_examine = 1;

% if constant choice is 1 specify the constant pH
constant_pH = 7;
% if constant choice is 2, specify the constant component molar value
constant_component_value = 0.02;

% Specify data to load. Should be the same as your gel_mechanics_model.m description
description = 'hNi';


%% Begin Program...
% Load data
cd(description)
load('structures.mat')

if y_choice == 1
	if constant_choice == 1
		% This is speciation at constant pH
		% first get index for constant pH. Just closest index right now.
		% unpack pH and find best index
		pH_list = speciation.pH_list;
		pH_diff_list = abs(constant_pH - pH_list);
		pH_match_list = find(pH_diff_list == min(pH_diff_list));
		pH_index = pH_match_list(1);
		
		plot_str = num2str(constant_pH);
		
		% get speciation data at this index. 
		for titration_number = 1:species_input_model.input.num_increments
			species(titration_number,:) = speciation.titration_number(titration_number).species(pH_index,:);
		end	
				
		species_names = [species_input_model.spec_names, species_input_model.spec_names_solids];
		model = species_input_model.model_all;
		
		%Prune out H column and OH column
		species(:,species_input_model.H_OH_indices) = [];
		species_names(species_input_model.H_OH_indices) = [];
		model(:,species_input_model.H_OH_indices) = [];
		
		% get number of ligands and number of metals
		% num_components, removing H
		num_components = species_input_model.number.components - 1;
		num_ligands = species_input_model.number.ligands;
		num_metals = species_input_model.number.metals;
		ligand_sum_index = num_ligands + 1;
		metal_sum_index = num_metals + 1;
		% with H component
		ligand_comp_indices_H = [num_metals+2:num_components+1];
		metal_comp_indices_H = [2:num_metals+1];
		% without H component
		ligand_comp_indices = [num_metals+1:num_components];
		metal_comp_indices = [1:num_metals];
		
		% figure out sum of ligand
		total_ligand = median(sum(species.*model(ligand_comp_indices_H,:),2));
		% figure out sum of metal
		total_metal = median(sum(species.*model(metal_comp_indices_H,:),2));
		total_metal = sum(species.*model(metal_comp_indices_H,:),2);
		
		molar_list = speciation.molar.component_concentration(:,component_number_examine);
		
		% Plot raw species
		figure
		plot(molar_list, species)
		legend(species_names)
		
		xlabel('Molar Concentration of Component')
		ylabel('Molar Concentration of Species')
		legend(species_names)
		axis square
		box on
		savefig(['speciation_',plot_str,'_species','.fig'])
		saveas(gcf,['speciation_',plot_str,'_species','.eps'],'epsc')
		
		% Plot of ligand raw
		plot_species_l_raw = species.*model(ligand_comp_indices_H,:);
		l_raw_sum = sum(plot_species_l_raw,1);
		l_raw_indices = find(l_raw_sum > 0);
		
		figure
		plot(molar_list, plot_species_l_raw)
		xlabel('Molar Concentration of Component')
		ylabel('Molar Concentration of Ligand')
		legend(species_names(l_raw_indices))
		axis square
		box on
		savefig(['speciation_',plot_str,'_l_raw','.fig'])
		saveas(gcf,['speciation_',plot_str,'_l_raw','.eps'],'epsc')
			
		% Plot of metal raw
		plot_species_m_raw = species.*model(metal_comp_indices_H,:);
		m_raw_sum = sum(plot_species_m_raw,1);
		m_raw_indices = find(m_raw_sum > 0);
		
		figure
		plot(molar_list, plot_species_m_raw)
		xlabel('Molar Concentration of Component')
		ylabel('Molar Concentration of Metal')
		legend(species_names(m_raw_indices))
		axis square
		box on
		savefig(['speciation_',plot_str,'_m_raw','.fig'])
		saveas(gcf,['speciation_',plot_str,'_m_raw','.eps'],'epsc')
		
		% Plot of ligand fraction
		plot_species_l_frac = species.*model(ligand_comp_indices_H,:) ./ total_ligand;
		l_frac_sum = sum(plot_species_l_frac,1);
		l_frac_indices = find(l_frac_sum > 0);
		
		figure
		plot(molar_list, plot_species_l_frac(:,l_frac_indices))
		xlabel('Molar Concentration of Component')
		ylabel('Fraction of Ligand')
		legend(species_names(l_frac_indices))
		axis square
		box on
		savefig(['speciation_',plot_str,'_l_frac','.fig'])
		saveas(gcf,['speciation_',plot_str,'_l_frac','.eps'],'epsc')
		
		% Plot of metal fraction
		plot_species_m_frac = species.*model(metal_comp_indices_H,:) ./ total_metal;
		m_frac_sum = sum(plot_species_m_frac,1);
		m_frac_indices = find(m_frac_sum > 0);
		
		figure
		plot(molar_list, plot_species_m_frac)
		xlabel('Molar Concentration of Component')
		ylabel('Fraction of Metal')
		legend(species_names(m_frac_indices))
		axis square
		box on
		savefig(['speciation_',plot_str,'_m_frac','.fig'])
		saveas(gcf,['speciation_',plot_str,'_m_frac','.eps'],'epsc')
		
		
	end
	
end

if y_choice == 1
	if constant_choice == 2
		% This is speciation at constant titration variable
		% first get index for constant titration. Just closest index right now.
		% unpack molar and find best index
		molar_list = speciation.molar.component_concentration(:,component_number_examine);
		molar_diff_list = abs(constant_component_value - molar_list);
		molar_match_list = find(molar_diff_list == min(molar_diff_list));
		molar_index = molar_match_list(1);
		
		plot_str = num2str(constant_component_value);
		
		% get speciation data at this index. 
		species = speciation.titration_number(molar_index).species;
		species_names = [species_input_model.spec_names, species_input_model.spec_names_solids];
		model = species_input_model.model_all;
		
		%Prune out H column and OH column
		species(:,species_input_model.H_OH_indices) = [];
		species_names(species_input_model.H_OH_indices) = [];
		model(:,species_input_model.H_OH_indices) = [];
		
		% get number of ligands and number of metals
		% num_components, removing H
		num_components = species_input_model.number.components - 1;
		num_ligands = species_input_model.number.ligands;
		num_metals = species_input_model.number.metals;
		ligand_sum_index = num_ligands + 1;
		metal_sum_index = num_metals + 1;
		% with H component
		ligand_comp_indices_H = [num_metals+2:num_components+1];
		metal_comp_indices_H = [2:num_metals+1];
		% without H component
		ligand_comp_indices = [num_metals+1:num_components];
		metal_comp_indices = [1:num_metals];
		
		% figure out sum of ligand
		total_ligand = median(sum(species.*model(ligand_comp_indices_H,:),2));
		% figure out sum of metal
		total_metal = median(sum(species.*model(metal_comp_indices_H,:),2));
		
		% get pH_list
		pH_list = speciation.pH_list;
		
		% Plot raw species
		figure
		plot(pH_list, species)
		legend(species_names)
		
		xlabel('pH')
		ylabel('Molar Concentration of Species')
		legend(species_names)
		axis square
		box on
		savefig(['speciation_',plot_str,'_species','.fig'])
		saveas(gcf,['speciation_',plot_str,'_species','.eps'],'epsc')
				
		% Plot of ligand raw
		plot_species_l_raw = species.*model(ligand_comp_indices_H,:);
		l_raw_sum = sum(plot_species_l_raw,1);
		l_raw_indices = find(l_raw_sum > 0);
		
		figure
		plot(pH_list, plot_species_l_raw)
		xlabel('pH')
		ylabel('Molar Concentration of Ligand')
		legend(species_names(l_raw_indices))
		axis square
		box on
		savefig(['speciation_',plot_str,'_l_raw','.fig'])
		saveas(gcf,['speciation_',plot_str,'_l_raw','.eps'],'epsc')
		
		% Plot of metal raw
		plot_species_m_raw = species.*model(metal_comp_indices_H,:);
		m_raw_sum = sum(plot_species_m_raw,1);
		m_raw_indices = find(m_raw_sum > 0);
		
		figure
		plot(pH_list, plot_species_m_raw)
		xlabel('pH')
		ylabel('Molar Concentration of Metal')
		legend(species_names(m_raw_indices))
		axis square
		box on
		savefig(['speciation_',plot_str,'_m_raw','.fig'])
		saveas(gcf,['speciation_',plot_str,'_m_raw','.eps'],'epsc')
		
		
		% Plot of ligand fraction
		plot_species_l_frac = species.*model(ligand_comp_indices_H,:) ./ total_ligand;
		l_frac_sum = sum(plot_species_l_frac,1);
		l_frac_indices = find(l_frac_sum > 0);
		
		figure
		plot(pH_list, plot_species_l_frac(:,l_frac_indices))
		xlabel('pH')
		ylabel('Fraction of Ligand')
		legend(species_names(l_frac_indices))
		axis square
		box on
		savefig(['speciation_',plot_str,'_l_frac','.fig'])
		saveas(gcf,['speciation_',plot_str,'_l_frac','.eps'],'epsc')
		
		% Plot of metal fraction
		plot_species_m_frac = species.*model(metal_comp_indices_H,:) ./ total_metal;
		m_frac_sum = sum(plot_species_m_frac,1);
		m_frac_indices = find(m_frac_sum > 0);
		
		figure
		plot(pH_list, plot_species_m_frac)
		xlabel('pH')
		ylabel('Fraction of Metal')
		legend(species_names(m_frac_indices))
		axis square
		box on
		savefig(['speciation_',plot_str,'_m_frac','.fig'])
		saveas(gcf,['speciation_',plot_str,'_m_frac','.eps'],'epsc')
		
	end
	
end


% Gp plots
if y_choice == 2
	if constant_choice == 1
			
		pH_list = speciation.pH_list;
		pH_diff_list = abs(constant_pH - pH_list);
		pH_match_list = find(pH_diff_list == min(pH_diff_list));
		pH_index = pH_match_list(1);
		
		molar_list = speciation.molar.component_concentration(:,component_number_examine);
		
		gp_matrix = mechanics.plot_gp.ligand(end).metal(end).Z;
		
		plot_str = num2str(constant_pH);
		
		figure
		plot(molar_list, gp_matrix(pH_index,:)./1000)
		xlabel('Molar Concentration of Component')
		ylabel('Plateau Modulus (kPa)')
		axis square
		box on
		savefig(['modulus_',plot_str,'_m_frac','.fig'])
		saveas(gcf,['modulus_',plot_str,'_m_frac','.eps'],'epsc')
		
		
	end
end


% Gp plots
if y_choice == 2
	if constant_choice == 2
		molar_list = speciation.molar.component_concentration(:,component_number_examine);
		molar_diff_list = abs(constant_component_value - molar_list);
		molar_match_list = find(molar_diff_list == min(molar_diff_list));
		molar_index = molar_match_list(1);
		
		pH_list = speciation.pH_list;
		
		gp_matrix = mechanics.plot_gp.ligand(end).metal(end).Z;
		
		plot_str = num2str(constant_component_value);
		
		figure
		plot(pH_list, gp_matrix(:,molar_index)./1000)
		xlabel('pH')
		ylabel('Plateau Modulus (kPa)')
		axis square
		box on
		savefig(['modulus_',plot_str,'_m_frac','.fig'])
		saveas(gcf,['modulus_',plot_str,'_m_frac','.eps'],'epsc')
		
		
	end
end


cd ..




















