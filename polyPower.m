function [output_poly] = polyPower(input_poly, power)

output_poly = 1;
	for k1 = 1:power
		output_poly = conv(output_poly, input_poly);
	end
end

