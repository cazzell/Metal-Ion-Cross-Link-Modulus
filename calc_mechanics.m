
function [mechanics] = calc_mechanics(species_input_model,speciation)

fprintf('calculating mechanics \n');

% Speciation output is in molar.

% Define some numbers
[~, num_titrations] = size(speciation.titration_number);
% num_components, removing H
num_components = species_input_model.number.components - 1;
num_ligands = species_input_model.number.ligands;
num_metals = species_input_model.number.metals;
num_pH = species_input_model.num_pH;

ligand_sum_index = num_ligands + 1;
metal_sum_index = num_metals + 1;

% with H component
ligand_comp_indices_H = [num_metals+2:num_components+1];
metal_comp_indices_H = [2:num_metals+1];

% without H component
ligand_comp_indices = [num_metals+1:num_components];
metal_comp_indices = [1:num_metals];

%% First thing to do next is to factor out the molar quantities that are in MOH complexes. MOH indices are not in species_input_model.

%% Calculate pl
% pm = 1 always
% Create pA/pB matrix
% Calculate fraction of ligand i and metal i reacted for each titration
% number and at each pH. Include a sum term by adding an extra ligand and
% metal indicating sum overall all ligands or metals respectively

for titration_number = 1:num_titrations
	% total molar cocentration of each ligand and metal and totals for each titration
	% these are the denominators for pA and pB
	mechanics.titration_number(titration_number).total_ligand = speciation.molar.component_concentration(titration_number, ligand_comp_indices);
	mechanics.titration_number(titration_number).total_metal = speciation.molar.component_concentration(titration_number, metal_comp_indices);
	
	mechanics.titration_number(titration_number).total_ligand(num_ligands+1) = sum(speciation.molar.component_concentration(titration_number, ligand_comp_indices));
	mechanics.titration_number(titration_number).total_metal(num_metals+1) = sum(speciation.molar.component_concentration(titration_number, metal_comp_indices));
	
	for ligand_number = 1:num_ligands+1
		if ligand_number > num_ligands
			ligand_row = sum(species_input_model.Model(ligand_comp_indices_H,:),1);
		else
			ligand_row = species_input_model.Model(ligand_comp_indices_H(ligand_number),:);
		end
		for metal_number = 1:num_metals+1
			if metal_number > num_ligands
				metal_row = sum(species_input_model.model_all(metal_comp_indices_H,:),1);
			else
				metal_row = species_input_model.model_all(metal_comp_indices_H(metal_number),:);
			end
			
			bound_indices = species_input_model.indices.bound.indices{ligand_number, metal_number};
			
			ligand_stoich = ligand_row(bound_indices);
			metal_stoich = metal_row(bound_indices);
			
			pL_numerator = nansum(ligand_stoich .* speciation.titration_number(titration_number).species(:,bound_indices),2);
			%pM_numerator = nansum(metal_stoich .* speciation.titration_number(titration_number).species(:,bound_indices),2);
			
			mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).pL_numerator = pL_numerator;
			%mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).pM_numerator = pM_numerator;
			mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).pL = pL_numerator ./ mechanics.titration_number(titration_number).total_ligand(ligand_number);
			%mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).pM = pM_numerator ./ mechanics.titration_number(titration_number).total_metal(metal_number);
			
			% Calculate amount of metal in MOH complexes
			
			%% right now this is forgetting about the solid complexes......
			% It should just be 1 for pB. Decide to keep or modify or prune
			% later.
			% 			MOH_indices = species_input_model.indices.MOH.indices{metal_number};
			%
			% 			metal_stoich = metal_row(MOH_indices);
			% 			MOH_sum = nansum(metal_stoich .* speciation.titration_number(titration_number).species(:,MOH_indices),2);
			%
			% 			mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).MOH_sum = MOH_sum;
			% 			mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).pB_neglect_MOH = pM_numerator ./ (mechanics.titration_number(titration_number).total_metal(metal_number) - MOH_sum);
			%
			% 			% Calculate amount of free metal???
			
			
			
		end
	end
	
	
	% 	for ligand_number = 1:num_ligands+1
	% 		for metal_number = 1:num_metals+1
	% 			temp_value = zeros(num_pH, 1);
	% 			MOH_indices = species_input_model.indices.MOH.indices{metal_number}
	% 			%for MOH_indicies = species_input_model.indices.MOH.indices{metal_number}
	% 			temp_value = nansum(speciation.titration_number(titration_number).species(:,MOH_indices),2);
	% 			%temp_value = speciation.titration_number(titration_number).species(:,MOH_indicies) + temp_value;
	% 			%end
	% 			if titration_number == 26
	% 				test = 5
	% 			end
	% 			mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).pB_neglect_MOH = mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).p_numerator ./ (mechanics.titration_number(titration_number).total_ligand(metal_number) - temp_value);
	% 		end
	% 	end
	
end

%% Calculate lf
% There is a seperate matrix for each ligand, and one with the ligands added together.
% The columns represent the functionality (column 2 is functionality = 2)
% The rows represent individual titration iterations

% Initialize the matrix for each ligand and combination
for ligand_number = 1:species_input_model.number.ligands+1
	mechanics.ligand(ligand_number).lf_molar = zeros(num_titrations,species_input_model.max_functionality+1);
	mechanics.ligand(ligand_number).lf_fraction = zeros(num_titrations,species_input_model.max_functionality+1);
end

% Calculate the molar quantities
for ligand_number = 1:species_input_model.number.ligands
	current_functionality = species_input_model.ligands_functionality(ligand_number);
	for titration_number = 1:num_titrations
		current_molar = speciation.molar.component_concentration(titration_number,ligand_number+species_input_model.number.metals);
		functionality_column = current_functionality;
		mechanics.ligand(ligand_number).lf_molar(titration_number,functionality_column) = current_molar;
		% Add to combined lf
		mechanics.ligand(species_input_model.number.ligands+1).lf_molar(titration_number,functionality_column) = mechanics.ligand(species_input_model.number.ligands+1).lf_molar(titration_number,functionality_column) + current_molar;
	end
end

% Calculate fractional values
for ligand_number = 1:species_input_model.number.ligands+1
	total_molar_ligand = sum(mechanics.ligand(ligand_number).lf_molar,2);
	mechanics.ligand(ligand_number).lf_fraction = mechanics.ligand(ligand_number).lf_molar ./ total_molar_ligand;
end


%% Calculate mg

for titration_number = 1:num_titrations
	mechanics.titration_number(titration_number).bound_L_molar = zeros(num_pH,species_input_model.number.metals+1);
	mechanics.titration_number(titration_number).bound_M_molar = zeros(num_pH,species_input_model.number.metals+1);
	% Initialize binding structure for the all metals (num_metals+1)
	for bind_state = 1:species_input_model.number.binding_states
		for metal = 1:species_input_model.number.metals + 1
			mechanics.titration_number(titration_number).metal(metal).bind_state(bind_state).molar_ligand = zeros(num_pH,1);
			mechanics.titration_number(titration_number).metal(metal).bind_state(bind_state).molar_metal = zeros(num_pH,1);
			
			mechanics.titration_number(titration_number).bind_state(bind_state).mg_by_metal_and_total = zeros(num_pH,species_input_model.number.metals+1);
			%mechanics.titration_number(titration_number).metal(metal).summations = zeros(num_pH,matrix_order);
		end
	end
	
	% Organize the molar volumes of metal and ligand in each binding state by metal and total
	for bind_state = 1:species_input_model.number.binding_states
		for ligand = 1:species_input_model.number.ligands
			for metal = 1:species_input_model.number.metals
				total_metal_index = species_input_model.number.metals + 1;
				current_index = species_input_model.indices.binding_state(bind_state).indices{ligand,metal};
				if ~isempty(current_index)
					mechanics.titration_number(titration_number).metal(metal).bind_state(bind_state).molar_ligand = (speciation.titration_number(titration_number).species(:,current_index) .* bind_state) + mechanics.titration_number(titration_number).metal(metal).bind_state(bind_state).molar_ligand;
					mechanics.titration_number(titration_number).metal(metal).bind_state(bind_state).molar_metal = (speciation.titration_number(titration_number).species(:,current_index)) + mechanics.titration_number(titration_number).metal(metal).bind_state(bind_state).molar_metal;
					
					mechanics.titration_number(titration_number).bound_L_molar(:,metal) = (speciation.titration_number(titration_number).species(:,current_index) .* bind_state) + mechanics.titration_number(titration_number).bound_L_molar(:,metal);
					mechanics.titration_number(titration_number).bound_M_molar(:,metal) = (speciation.titration_number(titration_number).species(:,current_index)) + mechanics.titration_number(titration_number).bound_M_molar(:,metal);
					
					% Now for the total index
					mechanics.titration_number(titration_number).metal(total_metal_index).bind_state(bind_state).molar_ligand = (speciation.titration_number(titration_number).species(:,current_index) .* bind_state) + mechanics.titration_number(titration_number).metal(total_metal_index).bind_state(bind_state).molar_ligand;
					mechanics.titration_number(titration_number).metal(total_metal_index).bind_state(bind_state).molar_metal = (speciation.titration_number(titration_number).species(:,current_index)) + mechanics.titration_number(titration_number).metal(total_metal_index).bind_state(bind_state).molar_metal;
					
					mechanics.titration_number(titration_number).bound_L_molar(:,total_metal_index) = (speciation.titration_number(titration_number).species(:,current_index) .* bind_state) + mechanics.titration_number(titration_number).bound_L_molar(:,total_metal_index);
					mechanics.titration_number(titration_number).bound_M_molar(:,total_metal_index) = (speciation.titration_number(titration_number).species(:,current_index)) + mechanics.titration_number(titration_number).bound_M_molar(:,total_metal_index);
					
				end
			end
		end
	end
	
	% Alternative p_l calculation
	% 	mechanics.titration(titration_number).p_A_by_metal_and_total = speciation.titration(titration_number).bound_L_molar./speciation.titration_parameters.total_poly_L_M_conc(titration_number);
	
	% Calculate mg data
	for bind_state = 1:species_input_model.number.binding_states
		for metal = 1:species_input_model.number.metals
			mechanics.titration_number(titration_number).bind_state(bind_state).mg_by_metal_and_total(:,metal) = mechanics.titration_number(titration_number).bind_state(bind_state).mg_by_metal_and_total(:,metal) + (mechanics.titration_number(titration_number).metal(metal).bind_state(bind_state).molar_metal ./ mechanics.titration_number(titration_number).bound_M_molar(:,end));
			mechanics.titration_number(titration_number).bind_state(bind_state).mg_by_metal_and_total(:,end) = mechanics.titration_number(titration_number).bind_state(bind_state).mg_by_metal_and_total(:,end) + (mechanics.titration_number(titration_number).metal(metal).bind_state(bind_state).molar_metal ./ mechanics.titration_number(titration_number).bound_M_molar(:,end));
		end
	end
end


%% Build Equation to solve

max_functionality = species_input_model.max_functionality;
max_binding = species_input_model.max_binding_state;

eqn_order = (max_functionality - 1) .* (max_binding-1);

% Intialize eqn

init_eqn = zeros(num_pH,eqn_order+1);
identity_vector = ones(num_pH,1);

for titration_number = 1:num_titrations
	for ligand_number = 1:num_ligands+1
		for metal_number = 1:num_metals+1
			
			eqn = init_eqn;
			% Easy stuff
			pL = mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).pL;
			eqn(:,end) = 1 - pL;
			eqn(:,end-1) = -1 .* identity_vector;
			
			for binding_state = 1:max_binding
				for functionality = 1:max_functionality
					current_lf = mechanics.ligand(ligand_number).lf_fraction(titration_number,functionality);
					raised_lf = current_lf .^ (binding_state - 1);
					
					current_mg = mechanics.titration_number(titration_number).bind_state(binding_state).mg_by_metal_and_total(:,metal_number);
					current_mg_lf = raised_lf .* current_mg;
					
					current_pL = mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).pL;
					
					current_pL_mg_lf = current_pL .* current_mg_lf;
					
					current_order = (functionality - 1) .* (binding_state - 1);
					current_column = eqn_order - current_order + 1;
					
					eqn(:,current_column) = eqn(:,current_column) + current_pL_mg_lf;
					
				end
			end
			mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).eqn = eqn;
		end
	end
end

%% Solve Equation
for titration_number = 1:num_titrations
	for ligand_number = 1:num_ligands+1
		for metal_number = 1:num_metals+1
			current_eqn = mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).eqn;
			for i = 1:num_pH
				output = roots(current_eqn(i,:));
				real_out= real(output);
				x = real_out(end);
				
				if (x<1) && (x>0)
					mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).x(i,1) = x;
				else
					mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).x(i,1) = 1;
				end
				
			end
		end
	end
end


%% Calculate crosslinks
%% Calculate molar crosslink concentrations

%% Arising from the polymers
% this is how we calculate the crosslinks for the polymer see eqn 45
max_crosslink_order = max([max_functionality, max_binding]);

for titration_number = 1:num_titrations
	for ligand_number = 1:num_ligands+1
		for metal_number = 1:num_metals+1
			% Initialize crosslink data
			mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).molar_crosslinks = zeros(num_pH, max_crosslink_order);
			mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).F_A_in = zeros(num_pH, 1);
		end
	end
end

for titration_number = 1:num_titrations
	for ligand_number = 1:num_ligands+1
		for metal_number = 1:num_metals+1
			current_functionalities = species_input_model.functionalities{ligand,metal};
			for f = current_functionalities
				if f > 2
					for m = [3:f]
						current_poly_functionality_molar = mechanics.ligand(ligand_number).lf_molar(titration_number,f) ./ f;
						F_A_out = mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).x;
						b = nchoosek(f,m) .* (F_A_out.^(f-m)) .* (1-F_A_out).^m;
						% need to multiply by molar volume
						b_molar = b .* current_poly_functionality_molar;
						mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).molar_crosslinks(:,m) = b_molar + mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).molar_crosslinks(:,m);
					end
				end
			end
		end
	end
end

%% Arising from metals with binding state > 3
% First need to calculate P(FAin)
for titration_number = 1:num_titrations
	for ligand_number = 1:num_ligands+1
		for metal_number = 1:num_metals+1
			current_functionalities = species_input_model.functionalities{ligand,metal};
			for f = current_functionalities
				lf_i = mechanics.ligand(ligand_number).lf_fraction(titration_number,f);
				F_A_out = mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).x;
				current_FAin_sum = lf_i.*(F_A_out.^(f-1));
				mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).F_A_in = current_FAin_sum + mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).F_A_in;
			end
		end
	end
end

% using eqn b.8 from Scott Grindy's thesis to calculate	crosslink concentration of metals
for titration_number = 1:num_titrations
	for ligand_number = 1:num_ligands+1
		for metal_number = 1:num_metals+1
			for binding_state = 1:species_input_model.number.binding_states
				if binding_state > 2
					mg = mechanics.titration_number(titration_number).bind_state(binding_state).mg_by_metal_and_total(:,metal_number);
					bound_M_molar = mechanics.titration_number(titration_number).bound_M_molar(:,metal_number);
					bound_M_molar_binding_state = bound_M_molar .* mg;
					F_A_in = mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).F_A_in;
					molar_crosslinks = bound_M_molar_binding_state .* (1-F_A_in).^(binding_state);
					mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).molar_crosslinks(:,binding_state) = molar_crosslinks + mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).molar_crosslinks(:,binding_state);
				end
			end
		end
	end
end

%
% Calculate concentration of crosslinks and chains in mol/m^3 (multiply molar quanitity by 1000
for titration_number = 1:num_titrations
	for ligand_number = 1:num_ligands+1
		for metal_number = 1:num_metals+1
			mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).total_mol_per_m_cubed_chains = zeros(num_pH, 1);
		end
	end
end

for titration_number = 1:num_titrations
	for ligand_number = 1:num_ligands+1
		for metal_number = 1:num_metals+1
			% Sum molar crosslinks and multiply by 1000 to get crosslinkn concentration
			mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).total_mol_per_m_cubed_crosslinks = sum(mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).molar_crosslinks,2).*1000;
			
			%[~,max_crosslink_order] = size(speciation.titration(titration_number).metal(metal).molar_crosslinks);
			
			% Calculate chain concentration
			for crosslink_type = 1:max_crosslink_order
				if crosslink_type > 2
					chains_per_crosslink_type = crosslink_type ./ 2;
					current_chains_vector = 1000 .* chains_per_crosslink_type .* mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).molar_crosslinks(:,crosslink_type);
					mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).total_mol_per_m_cubed_chains = current_chains_vector + mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).total_mol_per_m_cubed_chains;
				end
			end
		end
	end
end

% Calculate phantom network modulus
for titration_number = 1:num_titrations
	for ligand_number = 1:num_ligands+1
		for metal_number = 1:num_metals+1
			R = 8.314;
			T = 25+273.15;
			RT = R * T;
			mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).gp_phantom = RT .* (mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).total_mol_per_m_cubed_chains - mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).total_mol_per_m_cubed_crosslinks);
		end
	end
end

% Plot contour
for ligand_number = 1:num_ligands+1
	for metal_number = 1:num_metals+1
		[X,Y,Z] = plot_gp_contour(species_input_model,speciation,mechanics,metal_number,ligand_number);
		mechanics.plot_gp.ligand(ligand_number).metal(metal_number).X = X;
		mechanics.plot_gp.ligand(ligand_number).metal(metal_number).Y = Y;
		mechanics.plot_gp.ligand(ligand_number).metal(metal_number).Z = Z;
	end
end


