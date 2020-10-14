%Exercise - LMI

clear all
clf

%a)
A1 = [1 0.1;0.05 0.99];
A2 = [1 0.1;0 0];

B = [0;0.0787];

sdpvar Q(2,2)   % Square matrices are by default symmetric in Yalmip
sdpvar L(1,2)

F1 = [Q (Q*A1'+L'*B'); (A1*Q+B*L) Q];
F2 = [Q (Q*A2'+L'*B'); (A2*Q+B*L) Q];

Cstr = [F1>=0;F2>=0];

%out = optimize(Cstr,[]);    % Objective in second input argument.  Here we're only trying to fina a feasible solution.
                            % (empty objective function)

ops = sdpsettings('solver','mosek');
out = optimize(Cstr,[],ops);
                            
Qq = double(Q);             % Extracting solution
Ll = double(L);

P = inv(Qq);
K = Ll*P;

e1 = eig(A1+B*K)            % Checking solution
e2 = eig(A2+B*K)

% abs(eig) < 1, we have stable closed loop for both systems (and for
% arbitrarily fast time variation between the two systems)

%b)
% Q takes the role of inv(P)
C1 = [Q (A1+B*K)'*Q;Q*(A1+B*K) Q];
C2 = [Q (A2+B*K)'*Q;Q*(A2+B*K) Q];
Consb = [C1>=0;C2>=0];
f1 = [1 0];
f2 = [0 1];
Consb = [Consb;[1 f1*Q;Q*f1' Q]>= 0;[1 f2*Q;Q*f2' Q]>=0;[25 K*Q;Q*K' Q]>=0];
%out = optimize(Consb,-geomean(Q))
out = optimize(Consb,-geomean(Q),ops)

Qs = double(Q);

figure(1)
I = eye(2);
H = [I;-I;K;-K];
h = [ones(4,1);5;5];
Ph = Polyhedron(H,h);
plot(Ph,'Color','yellow') %Polytope fulfilling state and input constraints
hold on

sdpvar x(2,1)
plot(x'*inv(Qs)*x<=1)

% A larger ellipsoid could be found if K was designed to maximize the size
% of the ellipsoid.
%
%c)
% Let us try this in part c), following the description in slide 32.

% Only one input. Ej = (0,1), Ej- = (1 0)

sdpvar Y(1,2)
%Ej = 0:
C1 = [Q A1*Q+B*Y; Q*A1'+Y'*B' Q];
C2 = [Q A2*Q+B*Y; Q*A2'+Y'*B' Q];

%Ej = 1:
C3 = [Q A1*Q+B*L; Q*A1'+L'*B' Q];
C4 = [Q A2*Q+B*L; Q*A2'+L'*B' Q];

Consc = [C1>=0;C2>=0;C3>=0;C4>=0];
Consc = [Consc;[1 f1*Q;Q*f1' Q]>= 0;[1 f2*Q;Q*f2' Q]>=0;[25 Y;Y' Q]>=0];
out = optimize(Consc,-geomean(Q),ops)

%Extracting solution
Qq = double(Q);
L1 = double(L);
Y1 = double(Y);
P = inv(Qq);

K = L1*P;
Hv = Y1*P;

figure(2)

%Figure 2 shows in yellow the polytope defined by the state constraints and
%the constraint -umax < Kx < umax
%The inner ellipsoid is the ellipsoid found in b), the outer ellipsoid is
%the ellipsoid found in c).  Note that the outer ellipsoid goes outside the
%yellow box, since this ellipsoid shows the stable region with the
%_saturated_ input

I = eye(2);
H = [I;-I;K;-K];
h = [ones(4,1);5;5];
Ph = Polyhedron(H,h);
plot(Ph,'Color','yellow') %Polytope fulfilling state and input constraints
hold on

plot(x'*inv(Qq)*x<=1)

plot(x'*inv(Qs)*x<=1)

% Plotting trajectories

nstep = 200;
umax = 5;
t = linspace(1,nstep,nstep)';
for k = 1:18  %Initial condition index
    
    %Starting at 18 different positions on the border of the ellipse.
    %Figure 2 shows Model 1 trajectories in blue and Model 2 trajectories
    %in black
    
    %Figure 3 shows the saturated input
    % Model 1 input in blue and Model 2 input in black
    % We see that the saturated input adheres to the input constraint, also
    % for the initial point outside of the yellow box in Figure 2
    alpha = tan(k*20*2*pi/360);
    x1 = sqrt(1/(P(1,1)+2*alpha*P(1,2)+alpha^2*P(2,2)));
    if ((k>=5)&&(k<14))  %Choose negative solution
        x1 = -x1;
    end
    x2 = alpha*x1;
    x0 = [x1;x2];
    
    xt = trajcalc(x0,A1,B,K,umax,nstep);
    
    figure(2)
    plot(xt(:,1),xt(:,2),'b')
    hold on
    
    figure(3)
    plot(t,xt(:,3),'b')
    hold on
    
    % Model 2
    alpha = tan(k*20*2*pi/360);
    x1 = sqrt(1/(P(1,1)+2*alpha*P(1,2)+alpha^2*P(2,2)));
        if ((k>=5)&&(k<14))  %Choose negative solution
        x1 = -x1;
    end
    x2 = alpha*x1;
    x0 = [x1;x2];
    
    xt = trajcalc(x0,A2,B,K,umax,nstep);
    
    figure(2)
    plot(xt(:,1),xt(:,2),'k')
    hold on
    
    figure(3)
    plot(t,xt(:,3),'k')
    hold on
    
end
    
    
    
