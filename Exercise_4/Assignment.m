% Assignment for colloquium 4.  Stability analysis for Piecewise linear
% system

clear all

A1 = [ -0.1 1;-5 -0.1];
A2 = [ -0.1 5; -1 -0.1];

A3 = A1;
A4 = A2;

E1 = [-1 1;-1 -1]; %E1*x > 0 in region 1
E3 = -E1;
E2 = [-1 1;1 1];
E4 = -E2;

%Finding F's
%F's for neighbouring regions must be zero along rays on the boundary
%between regions.  E.g., F1*[-1;1] = F2*[-1 1], etc.
%Assembling these relationships as [F1 F2 F3 F4]*cstr = zeros(rf,ccstr]
cstr = [-1 0 0 1;1 0 0 1;1 1 0 0;-1 1 0 0; 0 -1 1 0;0 -1 -1 0; 0 0 -1 -1;0 0 1 -1];

%As dim(cstr) = (8,4), and cstr has full column rank, the maximum value for
%rf is 4. We can therefore fins Fi's of dim (4,2). F = [F1 F2 F3 F4] should
%lie in the left nullspace of cstr.  We want Fi's of the largest possible
%dimension, since this gives the most freedom in designing the lyapunov
%function.

[u,s,v] = svd(cstr);

%Four non-zero singular values. F is therefore found from the 4 last u
%vectors:

F= u(:,5:8)';
F1 = F(:,1:2);
F2 = F(:,3:4);
F3 = F(:,5:6);
F4 = F(:,7:8);

%May confirm that Fi*x = Fj*x for rays on boundaries between regions:
%F1*[-1;1]-F2*[-1; 1]
%F2*[1;1]-F3*[1; 1]
%F3*[1;-1]-F4*[1; -1]
%F4*[-1;-1]-F1*[-1; -1]

%Note that these Fi's are different from the one specified by Johansson and Rantzer.
%This is because any not-linearly-dependent set of combinations of the last
%4 column vectors of u could also be used for calculating F.
%I2 = eye(2)
%T = [E1' I2;E2' I2; E3' I2;E4' I2]
%T now contains the Fi's specified by Johansson and Ranzer
%We need to solve F'*X=T
%The least square solution is X = inv(F*F')*F*T
%X = inv(F*F')*F*T
%May check
%F'*X-T
%Which should be zero (up to numerical accuracy)

%We now have the Ei's and Fi's, and may proceed with stability analysis.

sdpvar T(4,4)
sdpvar U1(2,2)
sdpvar U2(2,2)
sdpvar U3(2,2)
sdpvar U4(2,2)
sdpvar W1(2,2)
sdpvar W2(2,2)
sdpvar W3(2,2)
sdpvar W4(2,2)



P1 = F1'*T*F1;  %Ensures continuity of LF across region boundaries
P2 = F2'*T*F2;
P3 = F3'*T*F3;
P4 = F4'*T*F4;

constr = [U1(:)>0,U2(:)>0,U3(:)>0,U4(:)>0,W1(:)>0,W2(:)>0,W3(:)>0,W4(:)>0]; % non-negative elements in the Ui's and Wi's

alpha = 0.200;  %Must iterate on alpha.  alpha = 0.2001 gives infeasible problem

constr = [constr,A1'*P1+P1*A1+E1'*U1*E1+alpha*P1<0]; %Derivative of LF < -alpha*P1 in region 1
constr = [constr,A2'*P2+P2*A2+E2'*U2*E2+alpha*P2<0];
constr = [constr,A3'*P3+P3*A3+E3'*U3*E3+alpha*P3<0];
constr = [constr,A4'*P4+P4*A4+E4'*U4*E4+alpha*P4<0];

constr = [constr,P1-E1'*W1*E1>0]; %LF>0 in region 1
constr = [constr,P2-E2'*W2*E2>0];
constr = [constr,P3-E3'*W3*E3>0];
constr = [constr,P4-E4'*W4*E4>0];

constr = [constr,trace(T)==4]; %(Superfluous) normalizing constraint

sol = optimize(constr,[]) %Try to find feasible solution (empty optimization criterion) for given alpha


%May check solution
% Td = double(T)
% U1d = double(U1)
% U2d = double(U2)
% U3d = double(U3)
% U4d = double(U4)
% W1d = double(W1)
% W2d = double(W2)
% W3d = double(W3)
% W4d = double(W4)





