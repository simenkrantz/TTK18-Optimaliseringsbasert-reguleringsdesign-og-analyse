clear
clc
close all
Tend=600;
Ts=1;
% Ts=100e-4;
A=[1  0.01;0.01 1];
B1=[0.001 0;0 -0.004];
b1=[0.09;0.09];
P=eye(2);
gamma= 295; %Starts reporting numerical problems  for gamma = 300
alpha = 0;  %Use alpha>0 for faster control, may then have to reduce gamma 
sdpvar x(2,1,'full') %t
sdpvar z(2,1,'full')
[u1,c11]=polynomial(x,1);
[u2,c12]=polynomial(x(2:2),1);
Cont = [u1 u2];
[u0,c0]=polynomial(x,2);
[s1,c2]=polynomial([x;z],2);
[s2,c3]=polynomial(x,2);
g=-(x'*P'*x-gamma);
Bx = B1*x;


F1 = [(1-alpha)*((u0+1)*P)  ((u0+1)*A+(Bx+b1)*Cont)'*P;
    P*((u0+1)*A+(Bx+b1)*Cont)     P*(u0+1)];
F=[sos([x;z]'*F1*[x;z]-s1*g),sos(u0),sos(s1),sos(s2),sos([(2^2)*(u0+1)-s2*g Cont*x;Cont*x (u0+1)])];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sedumi solver
%[sol,m,B,residual]=solvesos(F,[],sdpsettings('sos.model',1,'sos.postprocess',1,'sos.numblkdg',1e-4,'solver','sedumi','sedumi.eps',1e-11,'sedumi.numtol',1e-10),[c0;c11;c12;c2;c3]);

params = [c0;c11;c12;c2;c3];  %The parameters that are optimized
%sdpt3 solver
%[sol,m,B,residual]=solvesos(F,[],sdpsettings('sos.model',1,'sos.postprocess',1,'sos.numblkdg',1e-4,'solver','sdpt3','sdpt3.gaptol',1e-11,'sdpt3.inftol',1e-11,'sdpt3.steptol',1e-10,'sdpt3.maxit',100),params);

%MOSEK solver
%setting options
ops = sdpsettings('sos.model',1,'sos.postprocess',1,'sos.numblkdg',1e-4,'solver','mosek','mosek.MSK_DPAR_INTPNT_CO_TOL_DFEAS',1e-10,'mosek.MSK_DPAR_INTPNT_CO_TOL_PFEAS',1e-10,'mosek.MSK_DPAR_INTPNT_CO_TOL_REL_GAP',1e-10,'sos.numblk',0);
[sol,m,B,residual]=solvesos(F,[],ops,params);

% In newer Yalmip syntax, 'optimize' is used instead of 'solvesos' 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[h,parvar,Q,v]=sosd(F);


cc11=double(c11);  %Newer Yalmip syntax uses 'value' instead of 'double'.
cc12=double(c12);
cc0=double(c0);

%%%%%%


% CCont= [cc11(1)+x1*cc11(2)+x2*cc11(3)    cc12(1)+x2*cc12(2)];
% u00(x1,x2)=cc0(1)+x1*cc0(2)+x2*cc0(3)+x1^2*cc0(4)+x1*x2*cc0(5)+x2^2*cc0(6);
% u11(x1,x2)=CCont*[x1;x2];  %input definition
%plot the closed loop response
X(:,1)=[-10;13.9]; %initial condition
initial_cost=X(:,1)'*P*X(:,1)
for ii=1:Tend/Ts
    CCont= [cc11(1)+X(1,ii)*cc11(2)+X(2,ii)*cc11(3)    cc12(1)+X(2,ii)*cc12(2)];
    u00(ii)=cc0(1)+X(1,ii)*cc0(2)+X(2,ii)*cc0(3)+X(1,ii)^2*cc0(4)+X(1,ii)*X(2,ii)*cc0(5)+X(2,ii)^2*cc0(6);
    u11(ii)=CCont*[X(1,ii);X(2,ii)];
    ui(ii)=u11(ii)/(1+u00(ii));
    X(:,ii+1)=A*X(:,ii)+(B1*X(:,ii)+b1)*ui(ii);
    cost(ii)=X(:,ii)'*P*X(:,ii);
end

subplot(3,1,1)
plot(0:Tend/Ts-1,X(1,1:Tend/Ts))
hold on
plot(0:Tend/Ts-1,X(2,1:Tend/Ts),'r')
title('states evaluation in time, x1: blue, x2: red')
subplot(3,1,2)
plot(0:Tend/Ts-1,ui(1:Tend/Ts))
title('input')
subplot(3,1,3)
plot(0:Tend/Ts-1,cost(1:Tend/Ts))
title('cost function')

