fprintf('calculating mechanics');

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
		current_molar = speciation.molar.component_concentration(titration_number,ligand_number);
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
% 	%% Calculate molar crosslink concentrations
%
% 	% Arising from polymers with functionality > 2
% 	for metal = 1:species_input_model.number.metals+1
% 		speciation.titration(titration_number).metal(metal).molar_crosslinks = zeros(num_pH_titrations, max_crosslink_order);
% 		speciation.titration(titration_number).metal(metal).F_A_in = zeros(num_pH_titrations, 1);
% 	end
%
% 	for metal = 1:species_input_model.number.metals+1
% 		% this is how we calculate the crosslinks for the polymer see eqn 45
% 		functionality_vector = functionality_M_matrix(titration_number,:);
% 		functionality_index = 0;
% 		for current_functionality_molar = functionality_vector
% 			if current_functionality_molar ~= 0
% 				f = functionality_index;
% 				if f > 2
% 					for m = [3:f]
% 						F_A_out = speciation.titration(titration_number).metal(metal).F_A_out;
% 						b = nchoosek(f,m) .* (F_A_out.^(f-m)) .* (1-F_A_out).^m;
% 						% need to multiply by molar volume
% 						b_molar = b .* current_functionality_molar;
% 						speciation.titration(titration_number).metal(metal).molar_crosslinks(:,m) = b_molar + speciation.titration(titration_number).metal(metal).molar_crosslinks(:,m);
% 					end
% 				end
% 			end
% 			functionality_index = functionality_index + 1;
% 		end
% 	end
%
%
% 	% Arising from metals with binding state > 3
% 	% First need to calculate P(FAin)
% 	for metal = 1:species_input_model.number.metals+1
% 		f = 0;
% 		for af_i = speciation.titration_parameters.fraction_pol_L_funct_fplusone(titration_number,:)
% 			if f > 1
% 				F_A_out = speciation.titration(titration_number).metal(metal).F_A_out;
% 				current_FAin_sum = af_i.*(F_A_out.^(f-1));
% 				speciation.titration(titration_number).metal(metal).F_A_in = current_FAin_sum + speciation.titration(titration_number).metal(metal).F_A_in;
% 			end
% 			f = f+1;
% 		end
% 	end
%
% 	% using eqn b.8 from Scott Grindy's thesis to calculate	crosslink concentration of metals
% 	for metal = 1:species_input_model.number.metals+1
% 		for binding_state = 1:species_input_model.number.binding_states
% 			if binding_state > 2
% 				molar_metal_in_state = speciation.titration(titration_number).bind_state(bind_state).metal(metal).molar_metal;
% 				F_A_in = speciation.titration(titration_number).metal(metal).F_A_in;
% 				molar_crosslinks = molar_metal_in_state .* (1-F_A_in).^(binding_state);
% 				m = binding_state;
% 				speciation.titration(titration_number).metal(metal).molar_crosslinks(:,m) = molar_crosslinks + speciation.titration(titration_number).metal(metal).molar_crosslinks(:,m);
% 			end
% 		end
% 	end
%
% 	% Calculate concentration of crosslinks and chains in mol/m^3 (multiply
% 	% molar quanitity by 1000
% 	for metal = 1:species_input_model.number.metals+1
% 		speciation.titration(titration_number).metal(metal).total_mol_m_cubed_chains = zeros(num_pH_titrations, 1);
% 	end
%
% 	for metal = 1:species_input_model.number.metals+1
% 		speciation.titration(titration_number).metal(metal).total_mol_m_cubed_crosslinks = sum(speciation.titration(titration_number).metal(metal).molar_crosslinks,2).*1000;
% 		[~,max_crosslink_order] = size(speciation.titration(titration_number).metal(metal).molar_crosslinks);
% 		for crosslink_type = 1:max_crosslink_order
% 			if crosslink_type > 2
% 				chains_per_crosslink_type = crosslink_type ./ 2;
% 				current_chains_vector = 1000 .* chains_per_crosslink_type .* speciation.titration(titration_number).metal(metal).molar_crosslinks(:,crosslink_type);
% 				speciation.titration(titration_number).metal(metal).total_mol_m_cubed_chains = current_chains_vector + speciation.titration(titration_number).metal(metal).total_mol_m_cubed_chains;
% 			end
% 		end
% 	end
%
% 	for metal = 1:species_input_model.number.metals+1
% 		R = 8.314;
% 		T = 25+273.15;
% 		RT = R * T;
% 		speciation.titration(titration_number).metal(metal).gp_phantom = RT .* (speciation.titration(titration_number).metal(metal).total_mol_m_cubed_chains - speciation.titration(titration_number).metal(metal).total_mol_m_cubed_crosslinks);
% 	end
%
%
% 	%% Save Data From Each Titration Loop
% 	speciation.titration(titration_number).p_A_by_metal_and_total = speciation.titration(titration_number).bound_L_molar./speciation.titration_parameters.total_poly_L_M_conc(titration_number);
% 	speciation.titration(titration_number).species = C;
% 	speciation.titration(titration_number).pH = -log10(C(:,species_input_model.H_OH_indices(1)));
% 	speciation.titration(titration_number).component_total_M_concentrations = component_concentration_matrix(1:end-1,titration_number);
% 	speciation.titration_parameters.proton_total_M_concentrations = proton_conc';
%





%%



% %% Needed model data from generate_models.m
% % Model = species_input_model.Model;
% % beta = species_input_model.beta;
% max_functionality = species_input_model.max_functionality;
% max_g = species_input_model.max_binding_state;
% max_crosslink_order = max([max_functionality, max_g]);
%
% %% Titration parameters
% % Raster pH
% num_pH_titrations = 281;
% pH = linspace(species_input_model.pH(1),species_input_model.pH(2),281);
%
%
% % proton_conc_M_start = species_input_model.NaOH_range(1);
% % proton_conc_M_final = species_input_model.NaOH_range(2);
% % proton_conc = linspace(proton_conc_M_start, proton_conc_M_final, num_pH_titrations);
%
% %% Raster molar concentrations of components
% num_increments = 50;
% initial_M = species_input_model.initial_M;
% final_M = species_input_model.final_M;
% num_components_rastered = length(initial_M);
%
% component_concentration_matrix = zeros(num_increments, num_components_rastered);
% for column = 1:num_components_rastered
% 	for row = 1:num_increments
% 		linspace_vector = linspace(initial_M(column),final_M(column),num_increments);
% 		 component_concentration_matrix(:,column) = linspace_vector;
% 	end
% end
%
% %% Create af Matrix.
% % There is a seperate matrix for each ligand, and one with the ligands added together.
% % The columns represent the functionality, ordered highest to lowest (max functionality to 0)
% % The rows represent individual titration interations
%
% % Initialize the matrix for each ligand and combination
% for ligand_number = 1:species_input_model.number.ligands+1
% 	speciation.ligand(ligand_number).af_molar = zeros(num_increments,species_input_model.max_functionality+1);
% 	speciation.ligand(ligand_number).af_fraction = zeros(num_increments,species_input_model.max_functionality+1);
% end
%
% % Calculate the molar quantities
% for ligand_number = 1:species_input_model.number.ligands
% 	current_functionality = species_input_model.ligands_functionality(ligand_number);
% 	for titration_number = 1:num_increments
% 		current_molar = component_concentration_matrix(titration_number,ligand_number);
% 		functionality_column = species_input_model.max_functionality + 1 - current_functionality;
% 		speciation.ligand(ligand_number).af_molar(titration_number,functionality_column) = current_molar;
% 		% Add to combined ligand af
% 		speciation.ligand(species_input_model.number.ligands+1).af_molar(titration_number,functionality_column) = speciation.ligand(species_input_model.number.ligands+1).af_molar(titration_number,functionality_column) + current_molar;
% 	end
% end
%
% % Calculate fractional values
% for ligand_number = 1:species_input_model.number.ligands+1
% 	total_molar_ligand = sum(speciation.ligand(ligand_number).af_molar,2);
% 	speciation.ligand(ligand_number).af_fraction = speciation.ligand(ligand_number).af_molar ./ total_molar_ligand;
% end
%
% %% Solving for inner functionality summation which is independent of pH and metal identity
%
% % modify fraction to be the proper sum over f
% % sum_f of af*x^f-1 means we need to lower the term of everything by one
% for ligand_number = 1:species_input_model.number.ligands+1
% 	speciation.ligand(ligand_number).af_sum = speciation.ligand(ligand_number).af_fraction(:,1:end-1);
% end
%
% % determine the max order of the polynomial
% matrix_order = (max_functionality - 1) + (max_g - 1) + 1;
% % make initializing matrix
% initialized_inner_sum_matrix = zeros(num_increments,matrix_order);
%
% % Raise each term to g-1 power
% % Need to create a new matrix for each binding state
% for ligand_number = 1:species_input_model.number.ligands+1
% 	for bind_state = 1:species_input_model.number.binding_states
% 		%multiply af_sum by g-1
% 		%place in proper column
%
% 		% initialize inner_summation with matrix of zeros of proper dimensions
% 		speciation.ligand(ligand_number).bind_state(bind_state).inner_summation = initialized_inner_sum_matrix;
%
% 		% get the value of the current_af_sum (sum of (a_f * x^(f-1))
% 		current_af_sum = speciation.ligand(ligand_number).af_sum;
% 		% calculate (sum of (a_f^(g-1) * x^(f-1))
% 		modified_af_sum = current_af_sum .^ (bind_state - 1);
%
% 		% now need to properly place a_f^(g-1) into columns to get x^((f-1)+(g-1)
% 		[~,columns] = size(modified_af_sum);
% 		for active_column = 1:columns
% 			current_functionality = columns-active_column + 1;
% 			current_polynomial_order = current_functionality + bind_state - 2;
% 			placement_column = matrix_order - current_polynomial_order;
% 			speciation.ligand(ligand_number).bind_state(bind_state).inner_summation(:,placement_column) = modified_af_sum(:,active_column);
% 		end
% 	end
% end
%
%
%
%
%
%
% %% Everything above is independent of species???
%
%
% %% Run titration loop (calculate species concentrations as a function of pH)
%
% %%
% for titration_number=1:num_increments
% 	%% Speciation Calculation
% 	completed = int8((titration_number/num_increments) .* 100);
% 	fprintf('speciation calculation is %i%% completed \r \r', completed);
%
% 	current_comp_conc = component_concentration_matrix(titration_number,:);
%
% 	% Adding an extra entry here for the proton
% 	c_0 = [current_comp_conc,0];
%
% 	c_tot = repmat(c_0,num_pH_titrations,1);
% 	c_tot(:,end) = proton_conc;
%
% 	n_species = length(beta);
%
% 	% Initialize C, matrix of concentration in molar, for each species, for each titration
% 	C = zeros(num_pH_titrations,n_species);
%
% 	% c_comp_guess = [1e-10 1e-10 1e-10]; % initial guess for (L,M,H) in molar
% 	% Adding proton as a component that was not previously included in num_component
% 	c_comp_guess = ones(1,num_components_rastered+1);
% 	c_comp_guess = .1 .* c_comp_guess; % initial guess for (L,M,H) in molar
%
% 	for i=1:num_pH_titrations
% 		[C(i,:)]=NewtonRaphson(Model,beta,c_tot(i,:),c_comp_guess,i);
% 		c_comp_guess=C(i,1:num_components_rastered+1);
% 	end
%
% 	% Everything above here is good
%
%
% 	%% Generate Data for Gp calculation from speciation calculation
%
% 	%% Intialize Items
%
%
% 	%% Organize molar values ligand and metal in each binding state
% 	% speciation.titration(titration_number).ligand(ligand_number).metal(metal_number).bind_state(bind_state).molar_ligand=??
% 	% speciation.titration(titration_number).ligand(ligand_number).metal(metal_number).bind_state(bind_state).molar_metal=??
% 	% speciation.titration(titration_number).ligand(ligand_number).metal(metal_number).bind_state(bind_state).fraction_ligand=??
% 	% speciation.titration(titration_number).ligand(ligand_number).metal(metal_number).bind_state(bind_state).fraction_metal=??
% 	% fraction_metal gives molar metal for baseline so that pb = 1
% 	% fraction_ligand allows for calculation of pA
% 	% fraction of metal can also give bg
%
% 		%% Initialize a bunch of things
% 	speciation.titration(titration_number).bound_L_molar = zeros(num_pH_titrations,species_input_model.number.metals+1);
% 	speciation.titration(titration_number).bound_M_molar = zeros(num_pH_titrations,species_input_model.number.metals+1);
%
% 	% Initialize binding structure for the all metals (num_metals+1)
% 	for bind_state = 1:species_input_model.number.binding_states
% 		for metal = 1:species_input_model.number.metals + 1
% 			speciation.titration(titration_number).bind_state(bind_state).metal(metal).molar_ligand = zeros(num_pH_titrations,1);
% 			speciation.titration(titration_number).bind_state(bind_state).metal(metal).molar_metal = zeros(num_pH_titrations,1);
% 			speciation.titration(titration_number).bind_state(bind_state).bg_by_metal_and_total = zeros(num_pH_titrations,species_input_model.number.metals+1);
% 			speciation.titration(titration_number).metal(metal).summations = zeros(num_pH_titrations,matrix_order);
% 		end
% 	end
%
% 	% Organize the molar volumes of metal and ligand in each binding state by metal and total
% 	for bind_state = 1:species_input_model.number.binding_states
% 		for ligand = 1:species_input_model.number.ligands
% 			for metal = 1:species_input_model.number.metals
% 				total_metal_index = species_input_model.number.metals + 1;
% 				current_index = species_input_model.indices.poly_L_bound(bind_state).indices{ligand,metal};
% 				if ~isempty(current_index)
% 					speciation.titration(titration_number).bind_state(bind_state).metal(metal).molar_ligand = (C(:,current_index).* bind_state) + speciation.titration(titration_number).bind_state(bind_state).metal(metal).molar_ligand;
% 					speciation.titration(titration_number).bind_state(bind_state).metal(metal).molar_metal = (C(:,current_index)) + speciation.titration(titration_number).bind_state(bind_state).metal(metal).molar_metal;
% 					speciation.titration(titration_number).bound_L_molar(:,metal) = (C(:,current_index) .* bind_state) + speciation.titration(titration_number).bound_L_molar(:,metal);
% 					speciation.titration(titration_number).bound_M_molar(:,metal) = (C(:,current_index)) + speciation.titration(titration_number).bound_M_molar(:,metal);
% 					% Now for the total index
% 					speciation.titration(titration_number).bind_state(bind_state).metal(total_metal_index).molar_ligand = (C(:,current_index).* bind_state) + speciation.titration(titration_number).bind_state(bind_state).metal(total_metal_index).molar_ligand;
% 					speciation.titration(titration_number).bind_state(bind_state).metal(total_metal_index).molar_metal = (C(:,current_index)) + speciation.titration(titration_number).bind_state(bind_state).metal(total_metal_index).molar_metal;
% 				end
% 			end
% 		end
% 	end
% 	% Save this set of data
% 	speciation.titration(titration_number).bound_L_molar(:,species_input_model.number.metals+1) = sum(speciation.titration(titration_number).bound_L_molar(:,1:end-1),2);
% 	speciation.titration(titration_number).bound_M_molar(:,species_input_model.number.metals+1) = sum(speciation.titration(titration_number).bound_M_molar(:,1:end-1),2);
% 	speciation.titration(titration_number).p_A_by_metal_and_total = speciation.titration(titration_number).bound_L_molar./speciation.titration_parameters.total_poly_L_M_conc(titration_number);
%
% 	% Calculate bg data
% 	for bind_state = 1:species_input_model.number.binding_states
% 		for metal = 1:species_input_model.number.metals
% 			speciation.titration(titration_number).bind_state(bind_state).bg_by_metal_and_total(:,metal) = speciation.titration(titration_number).bind_state(bind_state).bg_by_metal_and_total(:,metal) + (speciation.titration(titration_number).bind_state(bind_state).metal(metal).molar_metal ./ speciation.titration(titration_number).bound_M_molar(:,end));
% 			speciation.titration(titration_number).bind_state(bind_state).bg_by_metal_and_total(:,end) = speciation.titration(titration_number).bind_state(bind_state).bg_by_metal_and_total(:,end) + (speciation.titration(titration_number).bind_state(bind_state).metal(metal).molar_metal ./ speciation.titration(titration_number).bound_M_molar(:,end));
% 		end
% 	end
%
%
% 	%% Generate F_A_out Matrix
% 	for bind_state = 1:species_input_model.number.binding_states
% 		for metal = 1:species_input_model.number.metals+1
% 			inner_summations = speciation.bind_state(bind_state).af_sum_full_matrix_to_gminusone(titration_number,:);
%
% 			bg = speciation.titration(titration_number).bind_state(bind_state).bg_by_metal_and_total(:,metal);
% 			speciation.titration(titration_number).bind_state(bind_state).metal(metal).all_summations = bg .* inner_summations;
%
% 			% sum over states
% 			speciation.titration(titration_number).metal(metal).summations = speciation.titration(titration_number).metal(metal).summations + speciation.titration(titration_number).bind_state(bind_state).metal(metal).all_summations;
% 		end
%
% 	end
%
% 	for metal = 1:species_input_model.number.metals+1
% 		pA = speciation.titration(titration_number).p_A_by_metal_and_total(:,metal);
% 		speciation.titration(titration_number).metal(metal).FA_out_eqn = speciation.titration(titration_number).metal(metal).summations .* pA;
% 		speciation.titration(titration_number).metal(metal).FA_out_eqn(:,end) = speciation.titration(titration_number).metal(metal).FA_out_eqn(:,end) - pA + 1;
% 		speciation.titration(titration_number).metal(metal).FA_out_eqn(:,end-1) = speciation.titration(titration_number).metal(metal).FA_out_eqn(:,end-1) - 1;
% 	end
%
%
%
% 	for metal = 1:species_input_model.number.metals+1
% 		polynomial_roots = zeros(num_pH_titrations,order_of_matrix-1);
% 		for i = 1:num_pH_titrations
% 			polynomial_roots(i,:) = roots(speciation.titration(titration_number).metal(metal).FA_out_eqn(i,:));
% 		end
% 		speciation.titration(titration_number).metal(metal).polynomial_roots = polynomial_roots;
% 		F_A_out = real(polynomial_roots(:,end));
% 		speciation.titration(titration_number).metal(metal).F_A_out = F_A_out;
% 	end
%
%
% 	%% Calculate molar crosslink concentrations
%
% 	% Arising from polymers with functionality > 2
% 	for metal = 1:species_input_model.number.metals+1
% 		speciation.titration(titration_number).metal(metal).molar_crosslinks = zeros(num_pH_titrations, max_crosslink_order);
% 		speciation.titration(titration_number).metal(metal).F_A_in = zeros(num_pH_titrations, 1);
% 	end
%
% 	for metal = 1:species_input_model.number.metals+1
% 		% this is how we calculate the crosslinks for the polymer see eqn 45
% 		functionality_vector = functionality_M_matrix(titration_number,:);
% 		functionality_index = 0;
% 		for current_functionality_molar = functionality_vector
% 			if current_functionality_molar ~= 0
% 				f = functionality_index;
% 				if f > 2
% 					for m = [3:f]
% 						F_A_out = speciation.titration(titration_number).metal(metal).F_A_out;
% 						b = nchoosek(f,m) .* (F_A_out.^(f-m)) .* (1-F_A_out).^m;
% 						% need to multiply by molar volume
% 						b_molar = b .* current_functionality_molar;
% 						speciation.titration(titration_number).metal(metal).molar_crosslinks(:,m) = b_molar + speciation.titration(titration_number).metal(metal).molar_crosslinks(:,m);
% 					end
% 				end
% 			end
% 			functionality_index = functionality_index + 1;
% 		end
% 	end
%
%
% 	% Arising from metals with binding state > 3
% 	% First need to calculate P(FAin)
% 	for metal = 1:species_input_model.number.metals+1
% 		f = 0;
% 		for af_i = speciation.titration_parameters.fraction_pol_L_funct_fplusone(titration_number,:)
% 			if f > 1
% 				F_A_out = speciation.titration(titration_number).metal(metal).F_A_out;
% 				current_FAin_sum = af_i.*(F_A_out.^(f-1));
% 				speciation.titration(titration_number).metal(metal).F_A_in = current_FAin_sum + speciation.titration(titration_number).metal(metal).F_A_in;
% 			end
% 			f = f+1;
% 		end
% 	end
%
% 	% using eqn b.8 from Scott Grindy's thesis to calculate	crosslink concentration of metals
% 	for metal = 1:species_input_model.number.metals+1
% 		for binding_state = 1:species_input_model.number.binding_states
% 			if binding_state > 2
% 				molar_metal_in_state = speciation.titration(titration_number).bind_state(bind_state).metal(metal).molar_metal;
% 				F_A_in = speciation.titration(titration_number).metal(metal).F_A_in;
% 				molar_crosslinks = molar_metal_in_state .* (1-F_A_in).^(binding_state);
% 				m = binding_state;
% 				speciation.titration(titration_number).metal(metal).molar_crosslinks(:,m) = molar_crosslinks + speciation.titration(titration_number).metal(metal).molar_crosslinks(:,m);
% 			end
% 		end
% 	end
%
% 	% Calculate concentration of crosslinks and chains in mol/m^3 (multiply
% 	% molar quanitity by 1000
% 	for metal = 1:species_input_model.number.metals+1
% 		speciation.titration(titration_number).metal(metal).total_mol_m_cubed_chains = zeros(num_pH_titrations, 1);
% 	end
%
% 	for metal = 1:species_input_model.number.metals+1
% 		speciation.titration(titration_number).metal(metal).total_mol_m_cubed_crosslinks = sum(speciation.titration(titration_number).metal(metal).molar_crosslinks,2).*1000;
% 		[~,max_crosslink_order] = size(speciation.titration(titration_number).metal(metal).molar_crosslinks);
% 		for crosslink_type = 1:max_crosslink_order
% 			if crosslink_type > 2
% 				chains_per_crosslink_type = crosslink_type ./ 2;
% 				current_chains_vector = 1000 .* chains_per_crosslink_type .* speciation.titration(titration_number).metal(metal).molar_crosslinks(:,crosslink_type);
% 				speciation.titration(titration_number).metal(metal).total_mol_m_cubed_chains = current_chains_vector + speciation.titration(titration_number).metal(metal).total_mol_m_cubed_chains;
% 			end
% 		end
% 	end
%
% 	for metal = 1:species_input_model.number.metals+1
% 		R = 8.314;
% 		T = 25+273.15;
% 		RT = R * T;
% 		speciation.titration(titration_number).metal(metal).gp_phantom = RT .* (speciation.titration(titration_number).metal(metal).total_mol_m_cubed_chains - speciation.titration(titration_number).metal(metal).total_mol_m_cubed_crosslinks);
% 	end
%
%
% 	%% Save Data From Each Titration Loop
% 	speciation.titration(titration_number).p_A_by_metal_and_total = speciation.titration(titration_number).bound_L_molar./speciation.titration_parameters.total_poly_L_M_conc(titration_number);
% 	speciation.titration(titration_number).species = C;
% 	speciation.titration(titration_number).pH = -log10(C(:,species_input_model.H_OH_indices(1)));
% 	speciation.titration(titration_number).component_total_M_concentrations = component_concentration_matrix(1:end-1,titration_number);
% 	speciation.titration_parameters.proton_total_M_concentrations = proton_conc';
%
% end
% end