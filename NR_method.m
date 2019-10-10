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

%% -------------------- NR method solids present
function [species,err,SI,solids]=NR_method(Asolution,Asolid,Ksolid,Ksolution,T,guess,iterations,criteria)

Nx=size(Asolution,2);
Ncp=size(Asolid,2);
Nc=size(Asolution,1);
X=guess;

for II=1:iterations
	
	
% 	if II < 100
% 		test = 30;
% 	end

	
	Xsolution=X(1:Nx);
	Xsolid=X(Nx+1:Nx+Ncp);
	
	logC=(Ksolution)+Asolution*log10(Xsolution);
	C=10.^(logC); % calc species
	
	Rmass=Asolution'*C+Asolid*Xsolid-T;
	
	Q=Asolid'*log10(Xsolution); SI=10.^(Q+Ksolid);
	RSI=ones(size(SI))-SI;
	
	% calc the jacobian
	z=zeros(Nx+Ncp,Nx+Ncp);
	
	for j=1:Nx
		for k=1:Nx
			for i=1:Nc
				z(j,k)=z(j,k)+Asolution(i,j)*Asolution(i,k)*C(i)/Xsolution(k);
			end
		end
	end
	for j=1:Nx
		for k=Nx+1:Nx+Ncp
			%t=Asolid';
			% z(j,k)=t(k-Nx,j);
			z(j,k)=Asolid(j,k-Nx);
		end
	end
	for j=Nx+1:Nx+Ncp
		for k=1:Nx
			z(j,k)=-1*Asolid(k,j-Nx)*(SI(j-Nx)/Xsolution(k));
		end
	end
	for j=Nx+1:Nx+Ncp
		for k=Nx+1:Nx+Ncp
			z(j,k)=0;
		end
	end
	
	R=[Rmass; RSI]; X=[Xsolution; Xsolid];
	
	deltaX=z\(-1*R);
	one_over_del=max([1, -1*deltaX'./(0.5*X')]); del=1/one_over_del;
	X=X+del*deltaX;
	
	tst=sum(abs(R));
	if tst<=criteria; break; end
end

logC=(Ksolution)+Asolution*log10(Xsolution);
C=10.^(logC); % calc species
RSI=ones(size(SI))-SI;

Rmass=Asolution'*C+Asolid*Xsolid-T;

err=[Rmass];
species=[C];
solids=Xsolid;

end