%% Code authored by Seth Allen Cazzell
% cazzell.lbi@gmail.com

% This is the program where you describe your system and input the
% variables.

% This is the master driver that calls other functions to generate
% theoretical predictions for metal-coordinated hydrogels of arbitrary
% functionality and combination of metals and ligands. Predicts gel plateau
% moduli as a function of pH and titrating concentrations.

% Please modify "generate_models.m" for your system of
% interest if desired ligand and metal is not present

%% Please reference
% "Expanding the stoichiometric window for metal cross-linked gel assembly using competition" PNAS, Cazzell (2019)
% "Engineering gelation in metal ion cross-linked hydrogels" MIT Thesis, Cazzell (2020)
% when appropriate.

%% Copyright 2020 Seth Allen Cazzell All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% This code calls functions modified or taken from the work of Donald Scott Smith,
% which is covered by the BSD 2 Clause License. Those functions include his BSD 2 license.

% Functions called in this code are directly or modified from Donald Scott Smith
% see his work
% Smith, D. Scott, "Solution of Simultaneous Chemical Equilibria in Heterogeneous Systems: Implementation in Matlab" (2019).
% Chemistry Faculty Publications. 14. https://scholars.wlu.ca/chem_faculty/14

% see reference 25 in "Expanding the stoichiometric window for metal cross-linked gel assembly using competition" PNAS, Cazzell (2019)
% cazzell.lbi@gmail.com

%% Setup and Clear Space
close all
clear
clc

% Add current path to working directory. Your current folder should be the
% main program folder with all of the .m files.
install_directory = pwd;
addpath(install_directory)

% Sometimes the matrix is initially unstable, results in a lot of outputs
% to the command window, which slows processing. These instabilities are
% corrected 
warning('off','all')

%% Input Variables
% description will be a directory that is created
% where everything will be saved
description = 'hNiCu';

ligands = {'h'};

% ligand functionality defines what the ligand is covalently attached to
% enter a value of
% 1 for a 1-arm polymer bound ligand or monofuntional small molecule
% 2 for a linear telechelic ligand
% f for a arbitrary f-arm star polymer ligand
ligands_functionality = [4];

metals = {'Ni','Cu'};

% Directly enforce M concentration of one ligand equivalence
lig_eqv = 0.04;

% Initial and final concentrations rastered across titration series
% General format is [M1, M2, Mn, L1, L2, Ln] or metals, then ligands
% Initial molar equivalence
initial_M_eqv = [1, 0, 1];
% Final molar equivalence
final_M_eqv =	[0, 1, 1];

% What is changing between intial and final molar specified above? ligand
% concentration? metal concentration? specify this here to serve as an
% x-axis label for the contour plots
label = 'Composition';

% Also, for plotting purposes, it's best if you can specify
% what the numerical range is for this variable
range = [0, 1];
% And specify the labels that you numerical labels that you want.
% You must have a label (an empty one is okay like {''}) at the beginning and end of
% the range
%x_tick_labels = [{''}, {'1/3'}, {'2/3'}, {'1'}, {'4/3'}, {''}];
x_tick_labels = [{'100% Ni'}, {'100% Cu'}];

% Fix the value of the maximum of the colorbar for the contour plot.
% This value is in kPa.
contour_max = 30;
% You can check if this is reasonable by checking the maximum after you
% plot with max(max(mechanics.plot_gp.ligand(end).metal(end).Z))./1000
% and then adjusting accordingly

% Resolution. Determine the resolution of the contour plot. The higher the
% resolution, the longer the calculations will take.
% y-axis, number of pH titrations
pH_range = [0,14];
num_pH = 251;
% x-axis, Number of titrations across initial to final concentrations
num_increments = 201;

% Use Hydroxides? Change this variable to 0 if you do not want to consider
% hydroxide competition. This is not advised, as hydroxide competition is
% present and more realistic. Value of 1 (default) includes hydroxide competition.
hydroxide_option = 1;


%% Package inputs
input.description = description;
input.ligands = ligands;
input.ligands_functionality = ligands_functionality;
input.metals = metals;
input.lig_eqv = lig_eqv;
input.initial_M_eqv = initial_M_eqv;
input.final_M_eqv =	final_M_eqv;
input.label = label;
input.range = range;
input.x_tick_labels = x_tick_labels;
input.pH_range = pH_range;
input.num_pH = num_pH;
input.num_increments = num_increments;
input.contour_max = contour_max;
input.hydroxide_option = hydroxide_option;

% Create and move into directory given by description
if ~exist(description, 'dir')
       mkdir(description)
end
cd(description)

%% Call Functions

% Generates model used to describe the desired system for the speciation program
[species_input_model] = generate_models(input);

% Runs speciation program to predict species concentrations vs. pH for each input titration
[speciation] = calc_speciation(species_input_model);

% Calculates and plots mechanics
[mechanics] = calc_mechanics(species_input_model,speciation);

% Code to save generated data
save('structures.mat','species_input_model','speciation','mechanics')

% Return to original directory
cd ..

