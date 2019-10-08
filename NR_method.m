% ?Solution of simultaneous chemical equilibria in heterogeneous systems: implementation in Matlab
% Scott Smith 2007
% Wilfrid Laurier University

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