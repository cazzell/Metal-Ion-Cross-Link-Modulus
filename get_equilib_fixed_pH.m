% ?Solution of simultaneous chemical equilibria in heterogeneous systems: implementation in Matlab
% Scott Smith 2007
% Wilfrid Laurier University

%% ----------- for fixed pH ----------------
function [Ksolution,Ksolid,Asolution,Asolid]=get_equilib_fixed_pH(KSOLUTION,KSOLID,ASOLUTION,ASOLID,pH)

% Transforms KSOLUTION, KSOLID, ASOLUTION and ASOLID for the current fixed pH value
% Removes the hydrogen component from A and rescales the formation constant

% Solutions
Ksolution=KSOLUTION-ASOLUTION(:,1)*pH;
% Strip away the proton column
[~,M]=size(ASOLUTION);
Asolution=[ASOLUTION(:,2:M)];

% Condensed Phases
Ksolid=KSOLID-ASOLID(:,1)*pH;
% Strip away the proton column
[~,M]=size(ASOLID);
Asolid=[ASOLID(:,2:M)];

end