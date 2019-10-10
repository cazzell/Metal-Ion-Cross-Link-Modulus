%% Copyright 2019 Donald Scott Smith All rights reserved. Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
% 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%% This code is taken directly from Donald Scott Smith
% see his work
% Smith, D. Scott, "Solution of Simultaneous Chemical Equilibria in Heterogeneous Systems: Implementation in Matlab" (2019).
% Chemistry Faculty Publications. 14. https://scholars.wlu.ca/chem_faculty/14

% see reference 25 in "Expanding the stoichiometric window for metal cross-linked gel assembly using competition" PNAS, Cazzell (2019)
% cazzell.lbi@gmail.com

%% ----------- NR method just solution species
function [species,err,SI]=NR_method_solution(Asolution,Asolid,Ksolid,Ksolution,T,guess,iterations,criteria)

Nx=size(Asolution,2);
Ncp=size(Asolid,1);
Nc=size(Asolution,1);
X=guess;

for II=1:iterations
	
	Xsolution=X(1:Nx);
	
	logC=(Ksolution)+Asolution*log10(Xsolution);
	C=10.^(logC); % calc species
	
	Rmass=Asolution'*C-T;
	
	Q=Asolid*log10(Xsolution); SI=10.^(Q+Ksolid);
	RSI=ones(size(SI))-SI;
	
	% calc the jacobian
	
	z=zeros(Nx,Nx);
	
	for j=1:Nx
		for k=1:Nx
			for i=1:Nc
				z(j,k)=z(j,k)+Asolution(i,j)*Asolution(i,k)*C(i)/Xsolution(k);
			end
		end
	end
	
	R=[Rmass]; X=[Xsolution];
	deltaX=z\(-1*R);
	one_over_del=max([1, -1*deltaX'./(0.5*X')]);
	del=1/one_over_del;
	X=X+del*deltaX;
	
	tst=sum(abs(R));
	if tst<=criteria; break; end
end
logC=(Ksolution)+Asolution*log10(Xsolution); C=10.^(logC); % calc species
RSI=ones(size(SI))-SI;

Rmass=Asolution'*C-T;
err=[Rmass];
species=[C];

end