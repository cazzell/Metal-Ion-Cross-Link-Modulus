%% Code authored by Seth Allen Cazzell
% cazzell.lbi@gmail.com

% This is the master driver that calls other functions to generate
% theoretical predictions for metal-coordinated hydrogels of arbitrary
% functionality and combination of metals and ligands. Predicts gel plateau
% moduli as a function of pH and titrating concentrations.

% This is the program that should be run to generate predictions. It calls
% the other functions. Please modify "generate_models.m" for your system of
% interest

%% Please reference
% "Expanding the stoichiometric window for metal cross-linked gel assembly using competition" PNAS, Cazzell (2019)
% when appropriate.

%% Copyright 2019 Seth Allen Cazzell All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
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

%%
close all
clear
clc

warning('off','all')

%% Call Functions

% Generates model used to describe the desired system for the speciation program
[species_input_model] = generate_models()

% Runs speciation program to predict species concentrations vs. pH for each input titration
[speciation] = calc_speciation(species_input_model)


% Calculates and plots mechanics
[mechanics] = calc_mechanics(species_input_model,speciation)

% Code to save generated data
save('structures.mat','species_input_model','speciation')

