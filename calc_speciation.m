%% Code authored by Seth Allen Cazzell
% cazzell.lbi@gmail.com

% This is the program to calculate the speciation of your system of interest.

%% Please reference
% "Expanding the stoichiometric window for metal cross-linked gel assembly using competition" PNAS, Cazzell (2019)
% when appropriate.

%% Copyright 2019 Seth Allen Cazzell All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% Copyright 2019 Donald Scott Smith All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% This code is modified from Donald Scott Smith
% see his work
% Smith, D. Scott, "Solution of Simultaneous Chemical Equilibria in Heterogeneous Systems: Implementation in Matlab" (2019).
% Chemistry Faculty Publications. 14. https://scholars.wlu.ca/chem_faculty/14

% see reference 25 in "Expanding the stoichiometric window for metal cross-linked gel assembly using competition" PNAS, Cazzell (2019)
% cazzell.lbi@gmail.com

%% Carrayrou et al AIChE journal 2002, 48, 894-904. % implement their method using their notation
% try HFO ppte as an example calc and at fixed pH
%function II=chem_equilib_Fe_fixedpH(pH,T)

%% For general method, see
% Practical Data Analysis in Chemistry
% Marcel Maeder & Yorck-Michael Neuhold
% University of Newcastle, Australia, 2006

%%
function [speciation] = calc_speciation(species_input_model)

%% Define Phase Space
% pH
num_pH_titrations = species_input_model.num_pH;
pH = linspace(species_input_model.pH(1), species_input_model.pH(2), num_pH_titrations);
speciation.pH_list = pH;

% Raster molar concentrations of components
num_increments = 100;
%num_increments = 10;
initial_M = species_input_model.initial_M;
final_M = species_input_model.final_M;
num_components_rastered = length(initial_M);

component_concentration_matrix = zeros(num_increments, num_components_rastered);
for column = 1:num_components_rastered
	for row = 1:num_increments
		linspace_vector = linspace(initial_M(column),final_M(column),num_increments);
		 component_concentration_matrix(:,column) = linspace_vector;
	end
end

speciation.molar.component_concentration = component_concentration_matrix;

speciation.SI.max_SI = zeros(1,size(species_input_model.calc_spec.ASOLID,1));
speciation.SI.bad_pixel_percentage = 0;
speciation.SI.bad_pixels = 0;

%% Run titration loop (calculate species concentrations as a function of pH)
for titration_number=1:num_increments
	%% Speciation Calculation
	completed = int8((titration_number/num_increments) .* 100);
	fprintf('speciation calculation is %i%% completed \r \r', completed);
	
	current_comp_conc = component_concentration_matrix(titration_number,:);
			
	T = current_comp_conc;
	
	% Extract Model Data from species_input_model
	KSOLUTION = species_input_model.calc_spec.KSOLUTION;
	KSOLID = species_input_model.calc_spec.KSOLID;
	ASOLUTION = species_input_model.calc_spec.ASOLUTION;
	ASOLID = species_input_model.calc_spec.ASOLID;
	SOLUTIONNAMES = species_input_model.calc_spec.SOLUTIONNAMES;
	SOLIDNAMES = species_input_model.calc_spec.SOLIDNAMES;
		
	% Make seed guess
	if titration_number > 1
		guess = speciation.titration_number(titration_number-1).species(1,2:num_components_rastered+1);
	else
		guess = ones(num_components_rastered,1)'.*(10^-10);
	end
		
	%numpt = number of points
	num_points=size(pH,2);
	%Ncp is number of condensed phases
	num_condensed_phases=size(ASOLID,1);
	% maximum number of times to run newton raphson
	iterations=100;
	% Acceptance criteria
	criteria=1e-16;
		
% 	% Initializes a variable named as the solid phases as a zero vector of length number of pH points
% 	for i=1:size(SOLIDNAMES,1)
% 		txt=[SOLIDNAMES(i,:),'=zeros(num_points,1);']; eval(txt)
% 	end
	
	% Loop at each pH value
	for pH_number=1:size(pH,2)
		
		% Recast equilibrium equations with fixed pH
		[Ksolution,Ksolid,Asolution,Asolid]=get_equilib_fixed_pH(KSOLUTION,KSOLID,ASOLUTION,ASOLID,pH(pH_number));
		
		Asolid_SI_check=Asolid;
		Ksolid_SI_check=Ksolid;
		
		% number of remaining components
		num_components = size(Asolution,2);
		% number of condensed phase species Ncp
		num_condensed_phases = size(Asolid,1);
		% number of solution species Nc
		num_solutions = size(Asolution,1);
		
		% calculate species using NR
		solids=zeros(1,num_condensed_phases);
		
						
		% Runs NR, but seeds guess from previous results when moving on to new pH values
		if pH_number==1
			[species,err,SI]=NR_method_solution(Asolution,Asolid,Ksolid,Ksolution,T',[guess(1:num_components)]',iterations,criteria);
		end
		if pH_number>1
			[species,err,SI]=NR_method_solution(Asolution,Asolid,Ksolid,Ksolution,T',[species(2:num_components+1)],iterations,criteria);
		end
				
		%% Runs calculation for condensed phase, makes corrections if saturation occurs
		for qq=1:num_condensed_phases
			
			[Y,I]=max(SI);
			
			if Y>1.000000001
				Iindex(qq)=I;
				% down selects Asolid and Ksolid to consider the worst
				% offender to deal with first. Once the program fixes this
				% problem, it will iterate again and choose the next worst
				% phase, adding it to Asolidtemp and Ksolidtemp to include
				% the previous progress. The program repeats this process
				% until all phases are added that have a saturation issue.
				Asolidtemp(qq,:)=Asolid_SI_check(I,:);
				Ksolidtemp(qq,:)=Ksolid_SI_check(I,:);
				
				% Just a guess to start things off for the first pH value
				solidguess(qq)=T(I)*0.5;
								
				% seeds with previous value
				if pH_number>1
					%txt=['solidguess(qq)=',SOLIDNAMES(I,:),'(pH_number-1);']; eval(txt);
					solidguess(qq) = solids(I);
				end
							
				guess=[species(2:num_components+1)' solidguess];
								
				% There was an error in the original program with the next
				% line of code. There was no index for the returning
				% solids. What I believe this code is meant to do, is first
				% try to fix the most problematic saturated species, then
				% add the next worst one etc. The problem with the original
				% code is that as Asolidtemp and Ksolidtemp added more and
				% more species, the returning container wasn't properly
				% sized to hold the solids that are returned.
						
				[species,err,SItst,solids(Iindex)]=NR_method(Asolution,Asolidtemp',Ksolidtemp,Ksolution,T',guess',iterations*1,criteria);
										
%				% Places values into solid named variable list		
% 				for q=1:size(solids,1)
% 					txt=[SOLIDNAMES(Iindex(q),:),'(pH_number)=solids(q);']; eval(txt)
% 				end
				
			end
			
			Q=Asolid*log10(species(2:num_components+1));
			SI=10.^(Q+Ksolid);
			Ifirst=I;
			
		end
		
		Q=Asolid*log10(species(2:num_components+1));
		SI=10.^(Q+Ksolid);
		SI_summary(pH_number,:)=SI;
		
		species_summary(pH_number,:)=species;
		mass_err_summary(pH_number,:)=(err(1));
		
		Asolidtemp=[]; Ksolidtemp=[];
		
		combined_species = [species',solids];
				
		speciation.titration_number(titration_number).species(pH_number,:) = combined_species;
		speciation.titration_number(titration_number).SI(pH_number,:) = SI;
		
		for SI_index = 1:length(SI)
			if SI(SI_index) > speciation.SI.max_SI(SI_index)
				speciation.SI.max_SI(SI_index) = SI(SI_index);
				speciation.SI.max_SI_titration_num(SI_index) = titration_number;
				speciation.SI.max_SI_pH_num(SI_index) = pH_number;
			end
			if SI(SI_index) > 1.01
				speciation.SI.bad_pixels = speciation.SI.bad_pixels + 1;
				speciation.SI.bad_pixel_percentage = 100 .* (speciation.SI.bad_pixels ./(num_condensed_phases .* num_pH_titrations .* num_increments));
			end
		end
		
	end

%	% Places values into solution species named variable list
% 	for i=1:size(species_summary,2)
% 		txt=[SOLUTIONNAMES(i,:),'=species_summary(:,i);']; eval(txt)
% 	end
	
end

end






