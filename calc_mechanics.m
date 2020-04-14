%% Code authored by Seth Allen Cazzell
% cazzell.lbi@gmail.com

% This is the program to calculate the mechanics from the speciation data

%% Please reference
% "Expanding the stoichiometric window for metal cross-linked gel assembly using competition" PNAS, Cazzell (2019)
% when appropriate.

%% Copyright 2019 Seth Allen Cazzell All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% For general method, see Miller and Macosko, and Grindy
% D. R. Miller, C. W. Macosko, A new derivation of post gel properties of network polymers. Macromolecules 9, 206-211 (1976).
% S. C. Grindy, "Complex mechanical design of bio-inspired model transient network hydrogels," PhD thesis, Massachusetts Institute of Technology (2017).

%%
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
end

%% Calculate lf
% There is a seperate matrix for each ligand, and one with the ligands added together.
% The columns represent the functionality (column 2 is functionality = 2)
% The rows represent individual titration iterations

% Initialize the matrix for each ligand and combination
for ligand_number = 1:num_ligands+1
	mechanics.ligand(ligand_number).lf_molar = zeros(num_titrations,species_input_model.max_functionality);
	mechanics.ligand(ligand_number).lf_fraction = zeros(num_titrations,species_input_model.max_functionality);
end

% Calculate the molar quantities
for ligand_number = 1:num_ligands
	current_functionality = species_input_model.ligands_functionality(ligand_number);
	for titration_number = 1:num_titrations
		current_molar = speciation.molar.component_concentration(titration_number,ligand_number+species_input_model.number.metals);
		functionality_column = current_functionality;
		mechanics.ligand(ligand_number).lf_molar(titration_number,functionality_column) = current_molar;
		% Add to combined lf
		mechanics.ligand(num_ligands+1).lf_molar(titration_number,functionality_column) = mechanics.ligand(num_ligands+1).lf_molar(titration_number,functionality_column) + current_molar;
	end
end

% Calculate fractional values
for ligand_number = 1:species_input_model.number.ligands+1
	total_molar_ligand = sum(mechanics.ligand(ligand_number).lf_molar,2);
	mechanics.ligand(ligand_number).lf_fraction = mechanics.ligand(ligand_number).lf_molar ./ total_molar_ligand;
end


%% Calculate mg

for titration_number = 1:num_titrations
	% Initialize binding structure for the all metals (num_metals+1)
	for bind_state = 1:species_input_model.number.binding_states
		for ligand_number = 1:species_input_model.number.ligands + 1
			for metal_number = 1:species_input_model.number.metals + 1
				mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bind_state(bind_state).molar_ligand = zeros(num_pH,1);
				mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bind_state(bind_state).molar_metal = zeros(num_pH,1);
				mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bind_state(bind_state).mg = zeros(num_pH,1);
			end
		end
	end
	% Intitialize bound molar and pL
	for ligand_number = 1:species_input_model.number.ligands + 1
		for metal_number = 1:species_input_model.number.metals + 1
			mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bound_M_molar = zeros(num_pH,1);
			mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bound_L_molar = zeros(num_pH,1);
			mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).pL = zeros(num_pH,1);
		end
	end
	
	% Organize the molar volumes of metal and ligand in each binding state by metal and total
	for bind_state = 1:species_input_model.number.binding_states
		for ligand_number = 1:num_ligands
			for metal_number = 1:num_metals
				total_metal_index = num_metals + 1;
				total_ligand_index = num_ligands + 1;
				current_index = species_input_model.indices.binding_state(bind_state).indices{ligand_number,metal_number};
				if ~isempty(current_index)

					mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bind_state(bind_state).molar_ligand = (speciation.titration_number(titration_number).species(:,current_index) .* bind_state) + mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bind_state(bind_state).molar_ligand;
					mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bind_state(bind_state).molar_metal = (speciation.titration_number(titration_number).species(:,current_index)) + mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bind_state(bind_state).molar_metal;
					
					mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bound_M_molar = (speciation.titration_number(titration_number).species(:,current_index)) + mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bound_M_molar;
					mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bound_L_molar = (speciation.titration_number(titration_number).species(:,current_index) .* bind_state) + mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bound_L_molar;
										
					% Total Metal
					mechanics.titration_number(titration_number).ligand(ligand_number).metal(total_metal_index).bound_M_molar = mechanics.titration_number(titration_number).ligand(ligand_number).metal(total_metal_index).bound_M_molar + (speciation.titration_number(titration_number).species(:,current_index));
					mechanics.titration_number(titration_number).ligand(ligand_number).metal(total_metal_index).bound_L_molar = mechanics.titration_number(titration_number).ligand(ligand_number).metal(total_metal_index).bound_L_molar + (speciation.titration_number(titration_number).species(:,current_index) .* bind_state);
					mechanics.titration_number(titration_number).ligand(ligand_number).metal(total_metal_index).bind_state(bind_state).molar_ligand = mechanics.titration_number(titration_number).ligand(ligand_number).metal(total_metal_index).bind_state(bind_state).molar_ligand + (speciation.titration_number(titration_number).species(:,current_index) .* bind_state);
					mechanics.titration_number(titration_number).ligand(ligand_number).metal(total_metal_index).bind_state(bind_state).molar_metal = mechanics.titration_number(titration_number).ligand(ligand_number).metal(total_metal_index).bind_state(bind_state).molar_metal + (speciation.titration_number(titration_number).species(:,current_index));
					% Total Ligand
					mechanics.titration_number(titration_number).ligand(total_ligand_index).metal(metal_number).bound_M_molar = mechanics.titration_number(titration_number).ligand(total_ligand_index).metal(metal_number).bound_M_molar + (speciation.titration_number(titration_number).species(:,current_index));
					mechanics.titration_number(titration_number).ligand(total_ligand_index).metal(metal_number).bound_L_molar = mechanics.titration_number(titration_number).ligand(total_ligand_index).metal(metal_number).bound_L_molar + (speciation.titration_number(titration_number).species(:,current_index) .* bind_state);
					mechanics.titration_number(titration_number).ligand(total_ligand_index).metal(metal_number).bind_state(bind_state).molar_ligand = mechanics.titration_number(titration_number).ligand(total_ligand_index).metal(metal_number).bind_state(bind_state).molar_ligand + (speciation.titration_number(titration_number).species(:,current_index) .* bind_state);
					mechanics.titration_number(titration_number).ligand(total_ligand_index).metal(metal_number).bind_state(bind_state).molar_metal = mechanics.titration_number(titration_number).ligand(total_ligand_index).metal(metal_number).bind_state(bind_state).molar_metal + (speciation.titration_number(titration_number).species(:,current_index));
					% Combined
					mechanics.titration_number(titration_number).ligand(total_ligand_index).metal(total_metal_index).bound_M_molar = mechanics.titration_number(titration_number).ligand(total_ligand_index).metal(total_metal_index).bound_M_molar + (speciation.titration_number(titration_number).species(:,current_index));
					mechanics.titration_number(titration_number).ligand(total_ligand_index).metal(total_metal_index).bound_L_molar = mechanics.titration_number(titration_number).ligand(total_ligand_index).metal(total_metal_index).bound_L_molar + (speciation.titration_number(titration_number).species(:,current_index)  .* bind_state);
					mechanics.titration_number(titration_number).ligand(total_ligand_index).metal(total_metal_index).bind_state(bind_state).molar_ligand = mechanics.titration_number(titration_number).ligand(total_ligand_index).metal(total_metal_index).bind_state(bind_state).molar_ligand + (speciation.titration_number(titration_number).species(:,current_index) .* bind_state);
					mechanics.titration_number(titration_number).ligand(total_ligand_index).metal(total_metal_index).bind_state(bind_state).molar_metal = mechanics.titration_number(titration_number).ligand(total_ligand_index).metal(total_metal_index).bind_state(bind_state).molar_metal + (speciation.titration_number(titration_number).species(:,current_index));
				end
			end
		end
	end
	
	for ligand_number = 1:num_ligands + 1
		for metal_number = 1:num_metals + 1
			mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).pL = mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bound_L_molar ./ mechanics.titration_number(titration_number).total_ligand(ligand_number);
		end
	end

	% Calculate mg data
	for bind_state = 1:species_input_model.number.binding_states
		for ligand_number = 1:num_ligands + 1
			for metal_number = 1:num_metals + 1
				mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bind_state(bind_state).mg = mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bind_state(bind_state).mg + (mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bind_state(bind_state).molar_metal ./ mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bound_M_molar);
			end
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
						
			build_eqn = init_eqn;
			f_eqn = init_eqn;
						
			for functionality = 1:max_functionality
				current_lf = mechanics.ligand(ligand_number).lf_fraction(titration_number,functionality);
				% This needs to be placed into f-1 position
				f_eqn(:,end-(functionality-1)) = f_eqn(:,end-(functionality-1)) + current_lf;
				% Here we would do all of the pM stuff if we decide we need
				% it. Currently ignoring because it is assumed to be 1.
			end
			
			for binding_state = 1:max_binding
				eqn = f_eqn;
				current_mg = mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bind_state(binding_state).mg;

				% Need to now raise the current equation to power g-1
				% This uses custom polyPower function to raise a vector
				% defining a polynomial to a power. The output vextor length 
				% can be a little weird, so have to do a little work on the
				% placement side to keep dimensions correct.
				for eqn_row = 1:length(eqn)
					out_poly = polyPower(eqn(eqn_row,:), binding_state - 1);
					for eqn_index = 1:eqn_order+1
						placement_index = eqn_index - 1;
						if eqn_index <= length(out_poly)
							eqn(eqn_row, end-placement_index) = out_poly(end-placement_index);
						end
					end
				end

				build_eqn = (eqn .* current_mg) + build_eqn;
				

			end
			
			current_pL = mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).pL;
			
			build_eqn = build_eqn .* current_pL;
				
			% Add additions that are outside the summations
			% 1 - pL term
			build_eqn(:,end) = build_eqn(:,end) + (1 - current_pL);
			% -x term
			build_eqn(:,end-1) = build_eqn(:,end-1) + (-1 .* identity_vector);
						
			mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).eqn = build_eqn;						
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
			current_functionalities = species_input_model.functionalities{ligand_number,metal_number};
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

%% Arising from metals with binding state > 2
% First need to calculate P(FAin)
for titration_number = 1:num_titrations
	for ligand_number = 1:num_ligands+1
		for metal_number = 1:num_metals+1
			current_functionalities = species_input_model.functionalities{ligand_number,metal_number};
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
					bound_M_molar_binding_state	= mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).bind_state(bind_state).molar_metal;
					F_A_in = mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).F_A_in;
					molar_crosslinks = bound_M_molar_binding_state .* (1-F_A_in).^(binding_state);
					mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).molar_crosslinks(:,binding_state) = molar_crosslinks + mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).molar_crosslinks(:,binding_state);
				end
			end
		end
	end
end

% Calculate concentration of crosslinks and chains in mol/m^3 (multiply molar quanitity by 1000
% Initialize
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
			% Sum molar crosslinks and multiply by 1000 to get crosslink concentration
			mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).total_mol_per_m_cubed_crosslinks = sum(mechanics.titration_number(titration_number).ligand(ligand_number).metal(metal_number).molar_crosslinks,2).*1000;
					
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
if num_ligands == 1 && num_metals == 1
	[X,Y,Z] = plot_gp_contour(species_input_model,speciation,mechanics,num_ligands+1,num_metals+1);
		mechanics.plot_gp.ligand(ligand_number).metal(metal_number).X = X;
		mechanics.plot_gp.ligand(ligand_number).metal(metal_number).Y = Y;
		mechanics.plot_gp.ligand(ligand_number).metal(metal_number).Z = Z;
else
	for ligand_number = num_ligands+1
		for metal_number = 1:num_metals+1
			[X,Y,Z] = plot_gp_contour(species_input_model,speciation,mechanics,metal_number,ligand_number);
			mechanics.plot_gp.ligand(ligand_number).metal(metal_number).X = X;
			mechanics.plot_gp.ligand(ligand_number).metal(metal_number).Y = Y;
			mechanics.plot_gp.ligand(ligand_number).metal(metal_number).Z = Z;
		end
	end
end

end


%% Nested Functions

%% polyPower
function [output_poly] = polyPower(input_poly, power)

output_poly = 1;
	for k1 = 1:power
		output_poly = conv(output_poly, input_poly);
	end
end

