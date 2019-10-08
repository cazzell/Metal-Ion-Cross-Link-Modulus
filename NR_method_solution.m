% ?Solution of simultaneous chemical equilibria in heterogeneous systems: implementation in Matlab
% Scott Smith 2007
% Wilfrid Laurier University

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