% Function to convert speciation data to plateau moduli predictions
% Based off of the methods developed by Miller and Macosko
% Miller, D. R. & Macosko, C. W. A New Derivation of Post Gel Properties of Network Polymers. Macromolecules 9, 206?211 (1976)
% Inspired by work of Scott Grindy
% Grindy, S. C. Complex Mechanical Design of Bio-Inspired Model Transient Network Hydrogels. (Massachusetts Institute of Technology, 2017)

function [speciation] = calc_Gp_linear_Miller_Macosko(speciation,species_input_model,metal_ratios,pH)

unbound_indices = species_input_model.unbound_L_indices;
mono_indices = species_input_model.mono_L_indices;
bis_indices = species_input_model.bis_L_indices;
tris_indices = species_input_model.tris_L_indices;

for metal_ratio_number = 1:length(metal_ratios)
    
    'MM'
    current_metal_ratio = metal_ratios(metal_ratio_number)
    %pH = speciation.metal_concentration(metal_ratio_number).pH;
    molar_ligand = speciation.molar_L;
    
    length_data = length(pH);
    
    fraction_empty = zeros(length_data,1);
    fraction_mono = zeros(length_data,1);
    fraction_bis = zeros(length_data,1);
    fraction_tris = zeros(length_data,1);
    
	
	
%     for indices = unbound_indices
%         fraction_empty = speciation.metal_concentration(metal_ratio_number).fraction.ligand(:,indices) + fraction_empty;
%     end
    
    for indices = mono_indices
        fraction_mono = speciation.metal_concentration(metal_ratio_number).fraction.ligand(:,indices) + fraction_mono;
    end
    
    for indices = bis_indices
        fraction_bis = speciation.metal_concentration(metal_ratio_number).fraction.ligand(:,indices) + fraction_bis;
    end
    
    for indices = tris_indices
        fraction_tris = speciation.metal_concentration(metal_ratio_number).fraction.ligand(:,indices) + fraction_tris;
    end
    
    %check_sum_L_frac = fraction_empty + fraction_mono + fraction_bis + fraction_tris;
    
	fraction_empty = 1 - fraction_mono - fraction_bis - fraction_tris;
	
    %initialize
    coeffs = zeros(length_data,3);
    
    %P(F_A_out)^2 term
    coeffs(:,1) = fraction_tris;
    %P(F_A_out)^1 term
    coeffs(:,2) = fraction_bis - 1;
    %P(F_A_out)^0 term
    coeffs(:,3) = fraction_mono + fraction_empty;
    
    polynomial_roots = zeros(length_data,2);
    
    for i = 1:length_data
        polynomial_roots(i,:) = roots(coeffs(i,:));
    end
    
    P_A_out = polynomial_roots(:,2);
    P_A_in = P_A_out;
    p_B_3 = (1/3) .* fraction_tris .* (1 - P_A_in).^3;
    
    % molar concentration of 3-functional crosslinks
    molar_X3 = p_B_3 .* molar_ligand;
    % convert to per m^3 from per L
    c_X3 = molar_X3 .* 1000;
    % c represents concentration in mol/m^3
    c_crosslinks = c_X3;
    % concentration of elastically active chains
    c_chains = 1.5*c_X3;
    
    R = 8.314;
    % Absolute Temperature
    T = 25+273.15;
    RT = R * T;
    %gp_affine = c_chains .* RT;
    gp_phantom = (c_chains - c_crosslinks) .* RT;
        
    speciation.metal_concentration(metal_ratio_number).gp_phantom = real(gp_phantom);
    
    %     speciation.metal_concentration(metal_ratio_number).(metal_ratio_number).pH = pH;
    %     speciation.metal_concentration(metal_ratio_number).(metal_ratio_number).gp_affine = gp_affine;
    %     speciation.metal_concentration(metal_ratio_number).(metal_ratio_number).fraction_empty = fraction_empty;
    %     speciation.metal_concentration(metal_ratio_number).(metal_ratio_number).fraction_mono = fraction_mono;
    %     speciation.metal_concentration(metal_ratio_number).(metal_ratio_number).fraction_bis = fraction_bis;
    %     speciation.metal_concentration(metal_ratio_number).(metal_ratio_number).fraction_tris = fraction_tris;
    %     speciation.metal_concentration(metal_ratio_number).(metal_ratio_number).check_sum = check_sum_L_frac;
        
end

end
