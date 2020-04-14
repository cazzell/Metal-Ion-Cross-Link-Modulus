%% Code authored by Seth Allen Cazzell
% cazzell.lbi@gmail.com

% This is the program to modify to generate predictions for your system of
% interest

%% Please reference
% "Expanding the stoichiometric window for metal cross-linked gel assembly using competition" PNAS, Cazzell (2019)
% when appropriate.

%% Copyright 2019 Seth Allen Cazzell All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%%
function [species_input_model] = generate_models(input)

% Unpack input
description = input.description;
ligands = input.ligands;
ligands_functionality = input.ligands_functionality;
metals = input.metals;
lig_eqv = input.lig_eqv;
initial_M_eqv = input.initial_M_eqv;
final_M_eqv = input.final_M_eqv;
pH_range = input.pH_range;
num_pH = input.num_pH;
hydroxide_option = input.hydroxide_option;

% Calculate initial and final molar concentration of components
% First make value of 0 = 1*10^-12, because 0 can break the calculation
initial_M_eqv(initial_M_eqv == 0) = 1*10^-6;
final_M_eqv(final_M_eqv == 0) = 1*10^-6;
initial_M = initial_M_eqv .* lig_eqv;
final_M = final_M_eqv .* lig_eqv;

% Toggles hydroxide competition
if hydroxide_option == 1
	hydroxide = {'OH'};
end

%% Formation Constants are taken from Smith, R. M. & Martell, A. E. "Critical Stability Constants," volume indicated

%% Protonation
% Stoichiometry
% Rows are Ligand
% Columns are Proton

% Nitrocatechol = n
% Volume 5
binding.n.H.matrix(1).constants(1,1) = 10.85;
binding.n.H.matrix(1).constants(1,2) = 6.69 + 10.85;

%Histamine = h
binding.h.H.matrix(1).constants(1,1) = 9.83;
binding.h.H.matrix(1).constants(1,2) = 6.07 + 9.83;

% Catechol = c
% Volume 5
binding.c.H.matrix(1).constants(1,1) = 12.8;
binding.c.H.matrix(1).constants(1,2) = 9.4 + 12.8;

% Imidazole = i
% Volume 6 page 249 (4-methylimidazole)
binding.i.H.matrix(1).constants(1,1) = 7.56;

% Hydroxide = OH
binding.OH.H.matrix(1).constants(1,1) = 13.997;

%% Metals-Ligand
% Stoichiometry
% Rows are Ligand
% Columns are Metal

% Nitrocatechol = n
% Iron = Fe
binding.n.Fe.matrix(1).constants(1,1) = 17.1;	% Volume 5
binding.n.Fe.matrix(1).constants(2,1) = 30.5;	% Volume 5
binding.n.Fe.matrix(1).constants(3,1) = 40;		% Volume 5
% Aluminum = Al
binding.n.Al.matrix(1).constants(1,1) = 13.74;	% Volume 6
binding.n.Al.matrix(1).constants(2,1) = 25.4;   % Volume 6
binding.n.Al.matrix(1).constants(3,1) = 34.3;   % Volume 6

% Histamine = h
% Nickel = Ni
binding.h.Ni.matrix(1).constants(1,1) = 6.78;	% Volume 2
binding.h.Ni.matrix(1).constants(2,1) = 11.78;	% Volume 2
binding.h.Ni.matrix(1).constants(3,1) = 14.9;	% Volume 2

% Cobolt = Co
binding.h.Co.matrix(1).constants(1,1) = 5.34;	% Volume 2
binding.h.Co.matrix(1).constants(2,1) = 9.09;	% Volume 2
binding.h.Co.matrix(1).constants(3,1) = 10.97;	% Volume 2

% Copper = Cu
binding.h.Cu.matrix(1).constants(1,1) = 9.56;	% Volume 2
binding.h.Cu.matrix(1).constants(2,1) = 16.13;	% Volume 2

% Zinc = Zn
binding.h.Zn.matrix(1).constants(1,1) = 5.03;	% Volume 2
binding.h.Zn.matrix(1).constants(2,1) = 9.81;	% Volume 2
binding.h.Zn.matrix(1).constants(3,1) = 12.09;	% Volume 2

% Catechol = c
% Iron = Fe
binding.c.Fe.matrix(1).constants(1,1) = 20.0;	% Volume 5
binding.c.Fe.matrix(1).constants(2,1) = 34.7;	% Volume 5
binding.c.Fe.matrix(1).constants(3,1) = 43.8;	% Volume 5

% Aluminum = Al
binding.c.Al.matrix(1).constants(1,1) = 16.3;	% Volume 3
binding.c.Al.matrix(1).constants(2,1) = 29.3;	% Volume 3
binding.c.Al.matrix(1).constants(3,1) = 37.6;	% Volume 3

% Imidazole = i (4-methylimidazole) Volume 6 page 249 and 250
% Ni
binding.i.Ni.matrix(1).constants(1,1) = 2.92;
binding.i.Ni.matrix(1).constants(2,1) = 5.25;
binding.i.Ni.matrix(1).constants(3,1) = 7.03;
binding.i.Ni.matrix(1).constants(4,1) = 8.25;

% Co
binding.i.Co.matrix(1).constants(1,1) = 2.34;
binding.i.Co.matrix(1).constants(2,1) = 4.09;
binding.i.Co.matrix(1).constants(3,1) = 5.33;
binding.i.Co.matrix(1).constants(4,1) = 6.67;

% Cu
binding.i.Cu.matrix(1).constants(1,1) = 4.18;
binding.i.Cu.matrix(1).constants(2,1) = 7.74;
binding.i.Cu.matrix(1).constants(3,1) = 10.70;
binding.i.Cu.matrix(1).constants(4,1) = 13.05;

% Zn
binding.i.Zn.matrix(1).constants(1,1) = 2.48;
binding.i.Zn.matrix(1).constants(2,1) = 5.06;
binding.i.Zn.matrix(1).constants(3,1) = 7.74;
binding.i.Zn.matrix(1).constants(4,1) = 10.52;


%% Metals-Hydroxide
% Stoichiometry
% Rows are Hydroxide
% Columns are Metal
% Multiple Matrices are if multiple species of the same stoichiometry form
% Last matrix is solids

% Hydroxide = OH
% Iron = Fe
binding.OH.Fe.matrix(1).constants(1,1) = -2.187;    % Volume 4
binding.OH.Fe.matrix(1).constants(2,1) = -5.694;    % Volume 4
binding.OH.Fe.matrix(2).constants(3,1) = -3.191;    % Volume 4
binding.OH.Fe.matrix(1).constants(4,1) = -21.588;	% Volume 4
binding.OH.Fe.matrix(1).constants(2,2) = -2.894;    % Volume 4
binding.OH.Fe.matrix(1).constants(4,3) = -6.288;    % Volume 4
% Solids
binding.OH.Fe.solid_matrix_number = 2;

% Aluminum = Al
binding.OH.Al.matrix(1).constants(1,1) = -4.987;	% Volume 4
binding.OH.Al.matrix(1).constants(2,1) = -9.294;	% Volume 4
binding.OH.Al.matrix(1).constants(3,1) = -14.991;	% Volume 4
binding.OH.Al.matrix(1).constants(4,1) = -22.998;	% Volume 4
binding.OH.Al.matrix(1).constants(2,2) = -7.694;	% Volume 4
binding.OH.Al.matrix(1).constants(4,3) = -13.888;	% Volume 4
binding.OH.Al.matrix(2).constants(3,1) = -8.491;	% Volume 4
% Solids
binding.OH.Al.solid_matrix_number = 2;

% Nickel = Ni
binding.OH.Ni.matrix(1).constants(1,1) = -9.897;	% Volume 4
binding.OH.Ni.matrix(1).constants(2,1) = -19.994;	% Volume 4
binding.OH.Ni.matrix(1).constants(3,1) = -30.991;	% Volume 4
binding.OH.Ni.matrix(1).constants(1,2) = -10.697;	% Volume 4
binding.OH.Ni.matrix(1).constants(4,4) = -27.688;	% Volume 4
binding.OH.Ni.matrix(2).constants(2,1) = -12.794;	% Volume 4
% Solids
binding.OH.Ni.solid_matrix_number = 2;

% Cobolt = Co
binding.OH.Co.matrix(1).constants(1,1) = -9.967;	% Volume 4
binding.OH.Co.matrix(1).constants(2,1) = -19.594;	% Volume 4
binding.OH.Co.matrix(1).constants(3,1) = -32.291;	% Volume 4
binding.OH.Co.matrix(1).constants(4,1) = -45.788;	% Volume 4
binding.OH.Co.matrix(1).constants(1,2) = -11.297;	% Volume 4
binding.OH.Co.matrix(1).constants(4,4) = -30.388;	% Volume 4
binding.OH.Co.matrix(2).constants(2,1) = -13.094;	% Volume 4
% Solids
binding.OH.Co.solid_matrix_number = 2;

% Copper = Cu
binding.OH.Cu.matrix(1).constants(1,1) = -7.697;	% Volume 4
binding.OH.Cu.matrix(1).constants(2,1) = -15.194;	% Volume 4
binding.OH.Cu.matrix(1).constants(3,1) = -27.491;	% Volume 4
binding.OH.Cu.matrix(1).constants(4,1) = -39.588;	% Volume 4
binding.OH.Cu.matrix(1).constants(2,2) = -10.294;	% Volume 4
binding.OH.Cu.matrix(1).constants(1,2) = -5.797;	% Volume 6
binding.OH.Cu.matrix(1).constants(4,3) = -20.788;	% Volume 6
binding.OH.Cu.matrix(2).constants(2,1) = -8.674;	% Volume 4
% Solids
binding.OH.Cu.solid_matrix_number = 2;

% Zinc = Zn
binding.OH.Zn.matrix(1).constants(1,1) = -8.997;	% Volume 4
binding.OH.Zn.matrix(1).constants(2,1) = -16.894;	% Volume 4
binding.OH.Zn.matrix(1).constants(3,1) = -28.391;	% Volume 4
binding.OH.Zn.matrix(1).constants(4,1) = -41.188;	% Volume 4
binding.OH.Zn.matrix(1).constants(1,2) = -8.997;	% Volume 4
binding.OH.Zn.matrix(1).constants(4,4) = -28.088;	% Volume 6
binding.OH.Zn.matrix(2).constants(2,1) = -12.474;	% Volume 4
% Solids
binding.OH.Zn.solid_matrix_number = 2;

%%
% DO NOT CHANGE WATER
proton = {'H'};
components = [(proton), (metals), (ligands)];
proton_index = 1;
component_numbers.hydroxide = length(components) + 1;
H_OH_indices = [proton_index, component_numbers.hydroxide];

% Calculate the number of each class of component
[~,number.components] = size(components);
[~,number.ligands] = size(ligands);
[~,number.metals] = size(metals);
[~,number.proton] = size(proton);

% Generate the row numbers for the components
component_numbers.proton = 1;
component_numbers.metals = 2:1+number.metals;
component_numbers.ligands = max(component_numbers.metals)+1:max(component_numbers.metals)+number.ligands;

base_vector = zeros(number.components,1);

% Make Identity Matrix for the Components
components_model = eye(number.components);
model = components_model;
model_log_beta = zeros(1,number.components);
% Add hydroxide for water self-ionization
hydroxide_base_vector = zeros(number.components,1);
hydroxide_base_vector(1) = -1;
model = [model, hydroxide_base_vector];
model_log_beta = [model_log_beta, -1.*binding.OH.H.matrix(1).constants(1,1)];

model_solids = [];
model_solids_log_beta = []; 

%% Get relevant constants for bound ligands
for ligand_number = 1:number.ligands
	%% Get Protonation
	[~,protonation_matrices] = size(binding.(ligands{ligand_number}).H.matrix);
	for protonation_matrix = 1:protonation_matrices
		vector = base_vector;
		input_matrix = binding.(ligands{ligand_number}).H.matrix(protonation_matrix).constants;
		[num_lig, num_protonation] = size(input_matrix);
		for row_num = 1:num_lig
			for col_num = 1:num_protonation
				constant = input_matrix(row_num,col_num);
				if constant ~= 0
					vector = base_vector;
					vector(component_numbers.ligands(ligand_number)) = row_num;
					vector(component_numbers.proton(1)) = col_num;
					log_beta = constant;
					model = [model,vector];
					model_log_beta = [model_log_beta,log_beta];
				end
			end
		end
	end
	
	%% Get Ligand-Metal Interactions
	for metal_number = 1:number.metals
		[~,metal_bound_ligand_matrices] = size(binding.(ligands{ligand_number}).(metals{metal_number}).matrix);
		for metal_matrix = 1:metal_bound_ligand_matrices
			vector = base_vector;
			input_matrix = binding.(ligands{ligand_number}).(metals{metal_number}).matrix(metal_bound_ligand_matrices).constants;
			[num_lig, num_metal] = size(input_matrix);
			for row_num = 1:num_lig
				for col_num = 1:num_metal
					constant = input_matrix(row_num,col_num);
					if constant ~= 0
						vector = base_vector;
						vector(component_numbers.ligands(ligand_number)) = row_num;
						vector(component_numbers.metals(metal_number)) = col_num;
						log_beta = constant;
						model = [model,vector];
						model_log_beta = [model_log_beta,log_beta];
					end
				end
			end
		end
	end
end

%% Get relevant constants for possible hydroxide complexes
% This can optionally be turned off by commenting out the hydroxide input
if exist('hydroxide','var')
	for metal_number = 1:number.metals
		[~,hydroxide_matrices] = size(binding.OH.(metals{metal_number}).matrix);
		for hydroxide_matrix = 1:hydroxide_matrices
			solid_matrix_number = binding.OH.(metals{metal_number}).solid_matrix_number;
			if hydroxide_matrix == solid_matrix_number
				vector = base_vector;
				input_matrix = binding.OH.(metals{metal_number}).matrix(hydroxide_matrix).constants;
				[num_OH, num_metals] = size(input_matrix);
				for row_num = 1:num_OH
					for col_num = 1:num_metals
						constant = input_matrix(row_num,col_num);
						if constant ~= 0
							vector = base_vector;
							vector(component_numbers.proton(1)) = -1 .* row_num;
							vector(component_numbers.metals(metal_number)) = col_num;
							log_beta = constant;
							model_solids = [model_solids,vector];
							model_solids_log_beta = [model_solids_log_beta,log_beta];
						end
					end
				end
			else
				vector = base_vector;
				input_matrix = binding.OH.(metals{metal_number}).matrix(hydroxide_matrix).constants;
				[num_OH, num_metals] = size(input_matrix);
				for row_num = 1:num_OH
					for col_num = 1:num_metals
						constant = input_matrix(row_num,col_num);
						if constant ~= 0
							vector = base_vector;
							vector(component_numbers.proton(1)) = -1 .* row_num;
							vector(component_numbers.metals(metal_number)) = col_num;
							log_beta = constant;
							model = [model,vector];
							model_log_beta = [model_log_beta,log_beta];
						end
					end
				end
			end
		end
	end
end

%%
% Get the indices for unbound, mono, bis and tris ligand.
[~,species_num] = size(model);

for ligand_number = 1:length(ligands)
	MOH_inc_all = 1;
	for metal_number = 1:length(metals)
		unbound_L_indices = [];
		bound_L_indices = [];
		MOH_indices = [];
		current_ligand_component_number = component_numbers.ligands(ligand_number);
		current_metal_component_number = component_numbers.metals(metal_number);
		current_ligand_row = model(current_ligand_component_number,:);
		current_metal_row = model(current_metal_component_number,:);
		H_row = model(1,:);
		unbound_inc = 1;
		bound_inc = 1;
		MOH_inc = 1;
		for current_species = 1:species_num
			L_value = current_ligand_row(current_species);
			M_value = current_metal_row(current_species);
			H_value = H_row(current_species);
			if M_value == 0
				if L_value > 0
					unbound_L_indices(unbound_inc) = current_species;
					unbound_inc = unbound_inc + 1;
				end
			end
			if M_value > 0
				if L_value > 0
					% Initialize species_input_model.indices.bind_state().indices for all ligands and metals
					if ~isfield(indices, 'binding_state') %checks if first time entering
						indices.binding_state(L_value).indices = cell(length(ligands), length(metals));
					elseif L_value > length(indices.binding_state)
						indices.binding_state(L_value).indices = cell(length(ligands), length(metals));
					end
					indices.binding_state(L_value).indices{ligand_number,metal_number} = current_species;
					bound_L_indices(bound_inc) = current_species;
					bound_inc = bound_inc + 1;
				end
				% Could add OH items here
				% Get MOH ligands. Only need to do once, so set ligand = 1
				if ligand_number == 1
					if H_value < 0
						MOH_indices(MOH_inc) = current_species;
						MOH_inc = MOH_inc + 1;
						MOH_indices_all(MOH_inc_all) = current_species;
						MOH_inc_all = MOH_inc_all + 1;
					end
				end
			end
			% save data to structure
			indices.unbound.indices{ligand_number,metal_number} = unbound_L_indices;
			functionalities{ligand_number,metal_number} = ligands_functionality(ligand_number);
			indices.bound.indices{ligand_number,metal_number} = bound_L_indices;
		end
		if ligand_number == 1
		%indices.MOH.indices{metal_number} = MOH_indices;
		%indices.MOH.indices{length(metals)+1} = MOH_indices_all;
		end
	end
end

long_hold_vector = [];
for ligand_number = 1:length(ligands)
	hold_vector = [];
	for metal_number = 1:length(metals)
		hold_vector = [hold_vector, indices.bound.indices{ligand_number,metal_number}];
		long_hold_vector = [long_hold_vector, indices.bound.indices{ligand_number,metal_number}];
	end
	indices.bound.indices{ligand_number,length(metals)+1} = hold_vector;
end
indices.bound.indices{length(ligands)+1,length(metals)+1} = long_hold_vector;

for metal_number = 1:length(metals)
	hold_vector = [];
	for ligand_number = 1:length(ligands)
		hold_vector = [hold_vector, indices.bound.indices{ligand_number,metal_number}];
	end
	indices.bound.indices{length(ligands)+1,metal_number} = hold_vector;
end

% Name all of the species
for current_species = 1:size(model,2)
	current_species_vector = model(:,current_species);
	index = 1;
	for current_component_number = 1:length(current_species_vector)
		current_stoich = current_species_vector(current_component_number);
		if current_stoich > 0
			component_name = components(current_component_number);
			component_stoich = num2str(current_stoich);
			if index == 1
				species_name = [component_name{1},'_', component_stoich];
			end
			if index > 1
				species_name = [species_name, '__',component_name{1},'_', component_stoich];
			end
			index = index + 1;
		end
		if current_stoich < 0
			component_name = components(current_component_number);
			component_stoich = num2str(-1.*current_stoich);
			if index == 1
				species_name = ['OH','_', component_stoich];
			end
			if index > 1
				species_name = [species_name, '__','OH','_', component_stoich];
			end
			index = index + 1;
		end
	end
	species_names{1,current_species} = species_name;
end

% Name solid species
for current_species = 1:size(model_solids,2)
	current_species_vector = model_solids(:,current_species);
	index = 1;
	for current_component_number = 1:length(current_species_vector)
		current_stoich = current_species_vector(current_component_number);
		if current_stoich > 0
			component_name = components(current_component_number);
			component_stoich = num2str(current_stoich);
			if index == 1
				species_name_solids = [component_name{1},'_', component_stoich];
			end
			if index > 1
				species_name_solids = [species_name_solids, '__',component_name{1},'_', component_stoich];
			end
			index = index + 1;
		end
		if current_stoich < 0
			component_name = components(current_component_number);
			component_stoich = num2str(-1.*current_stoich);
			if index == 1
				species_name_solids = ['OH','_', component_stoich];
			end
			if index > 1
				species_name_solids = [species_name_solids, '__','OH','_', component_stoich];
			end
			index = index + 1;
		end
	end
	species_names_solids{1,current_species} = species_name_solids;
end

% Add functionalities for combination
[rows,columns] = size(functionalities);

% Add to column
for row_num = 1:rows
	col_num = columns+1;
	functionalities{row_num,col_num} =  functionalities{row_num,col_num-1};
end

% Add to rows
for col_num = 1:columns+1
	functionality_vector = null(1,1);
	for row_num = 1:rows
		current_functionality = functionalities{row_num,col_num};
		if  ~ismember(current_functionality, functionality_vector)
			functionality_vector = [functionality_vector,current_functionality];
		end
	end
	functionalities{row_num+1,col_num} = functionality_vector;
end

species_input_model.description = description;
species_input_model.pH = pH_range;

species_input_model.spec_names = species_names;
species_input_model.Model = model;
species_input_model.log_beta = model_log_beta;

if exist('hydroxide','var')
	species_input_model.spec_names_solids = species_names_solids;
	species_input_model.Model_solids = model_solids;
	species_input_model.log_beta_solids = model_solids_log_beta;
	species_input_model.calc_spec.SOLIDNAMES = char(species_names_solids);
	species_input_model.calc_spec.KSOLID = model_solids_log_beta';
	species_input_model.calc_spec.ASOLID = model_solids';
else
	species_input_model.spec_names_solids = {};
	species_input_model.Model_solids = double.empty(number.components,0);
	species_input_model.log_beta_solids = double.empty(1,0);
	species_input_model.calc_spec.SOLIDNAMES = '';
	species_input_model.calc_spec.KSOLID = double.empty(0,1);
	species_input_model.calc_spec.ASOLID = double.empty(0,number.components);
end

species_input_model.calc_spec.KSOLUTION = model_log_beta';
species_input_model.calc_spec.ASOLUTION = model';
species_input_model.calc_spec.SOLUTIONNAMES = char(species_names);

species_input_model.components = components;
species_input_model.indices = indices;
species_input_model.functionalities = functionalities;
species_input_model.ligands_functionality = ligands_functionality;
species_input_model.max_functionality = max(ligands_functionality);
[~,species_input_model.max_binding_state] = size(indices.binding_state);
species_input_model.initial_M = initial_M;
species_input_model.final_M = final_M;
species_input_model.H_OH_indices = H_OH_indices;
species_input_model.number = number;
species_input_model.num_pH = num_pH;
species_input_model.ligands = ligands;
species_input_model.model_all = [model, model_solids];
species_input_model.log_beta_all = [model_log_beta, model_solids_log_beta];
% Saves other input variables for access later
species_input_model.input = input;

[~ , max_binding_state] = size(species_input_model.indices.binding_state);
species_input_model.number.binding_states = max_binding_state;

end
