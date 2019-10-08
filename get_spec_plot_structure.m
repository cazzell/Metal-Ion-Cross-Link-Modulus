% This is a function that returns a data structure organizing speciation data as a fraction of a given component.

function [spec_plot] = get_spec_plot_structure(row_index,speciation,species_input_model)

species_names = species_input_model.spec_names;
model = species_input_model.Model;
model_row = model(row_index,:);

% Get the species that contain the component you are interested in
number_of_species = 1;
for species_number = 1:length(model_row)
    if model_row(species_number) > 0.1
        species_index(number_of_species) = species_number;
        number_of_species = number_of_species + 1;
    end
end

species_names = species_names(species_index);

% For each experimentally relevant titration, organize speciation data
metal_plot_index = 1;
for current_titration_index = speciation.critical_titration.indices
    spec_plot.metal_conc(metal_plot_index).pH = speciation.pH;
    %spec_plot.metal_conc(metal_plot_index).v_added = speciation.v_added;
    %spec_plot.metal_conc(metal_plot_index).norm_v_added = speciation.norm_v_added;
    spec_plot.metal_conc(metal_plot_index).species_names = species_names;
	if row_index == 3 
		spec_plot.metal_conc(metal_plot_index).species = speciation.metal_concentration(current_titration_index).fraction.ligand(:,species_index);
	end
	if row_index == 2
		spec_plot.metal_conc(metal_plot_index).species = speciation.metal_concentration(current_titration_index).fraction.metal(:,species_index);
	end
    %spec_plot.metal_conc(metal_plot_index).species = speciation.metal_concentration(current_titration_index).fraction_row(row_index).frac(:,species_index);
    %spec_plot.metal_conc(metal_plot_index).check_sum = sum(spec_plot.metal_conc(metal_plot_index).species,2);
    metal_plot_index = metal_plot_index + 1;
end

end

