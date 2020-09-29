%Tentative solution to voluntary exercise
nx = 2; %number of states
nu = 1; %number of inputs

A = [0.7326 -0.0861;0.1722 0.9909];
B = [0.0609;0.0064];
ncx = 4; % number of state constraints
H = [eye(2);-eye(2)];
h = 100*ones(ncx,1);
ncu = 2; % number of input constraints
Hu = [1;-1];
hu = 2*ones(ncu,1); %Change H,h,Hu,hu when calculating contraction factor for a different set.  Ml and Mu (below) may also have to change.

%Variable organization:  v = [x; u; \gamma; \lambda; s] 
%(states, input,contraction factor, Lagrange multipliers, binary variables denoting active constraints)

Ml = diag([0.01 0.01 0.01 0.01 1 1]); % The first equality constraint ensures (in this case) that h'*\lam1 = 1, where \lam1 are the Lagrangians corresdponding to the contraction constraints,
                      %giving an upper bound to the Lagrangians
Mu = diag([300 300 300 300 6 6]); %For larger problems it could be an advantage to 
                        %find individual upper bounds on constraints by
                        %solving a series of LPs.
                        %Here I have simply adjusted manually until the
                        %corresponding constraints are no longer active at
                        %the optimum.  Remember: too small values will give
                        %you the wrong solution.  Too large values may give
                        %inaccurate solution due to ill-conditioning.
                        
big = 1000;

blx = [-100*[1;1];-2;0;zeros(ncx+ncu,1);zeros(ncx+ncu,1)];  %Lower bound on variables
bux = [100*[1;1]; 2; 2;ones(ncx+ncu,1);ones(ncx+ncu,1)]; %Upper bound on variables

Aineq = zeros(3*ncx+2*ncu,nx+nu+1+2*(ncx+ncu)); %Initializing matrix of inequality constraints
bineq = zeros(3*ncx+2*ncu,1);

Aeq = zeros(1+nu,nx+nu+1+2*(ncx+ncu)); %Initializing matrix of equality constraints
beq = zeros(1+nu,1);

np = 0;
Aineq(np+1:np+ncx,1:nx+nu+1) = [H*A H*B -h]; % H(Ax+Bu) \le \gamma*h
%RHS already initializer to zero

np = np + ncx;


Aineq(np+1:np+ncx+ncu,nx+nu+2:nx+nu+1+2*(ncx+ncu)) = [eye(ncx+ncu) -Ml]; % \lambda \le Ml*s
% RHS already zero

np = np+ncx+ncu;
Aineq(np+1:np+ncx+ncu,:) = [[-H*A -H*B h zeros(ncx,ncx+ncu);zeros(ncu,nx) -Hu zeros(ncu,1) zeros(ncu,ncx+ncu)] Mu];
bineq(np+1:np+ncx+ncu,1) = [zeros(ncx,1);-hu]+Mu*ones(ncx+ncu,1);      %[ H*A*x+H*B*u \ge h-Mu*(1-s); Hu*u \ge hu-Mu(1-s)]

Aeq(1,nx+nu+2:nx+nu+1+ncx) = h'; % Derivative of Lagrangian w.r.t. \gamma = 0
beq(1,1) = 1;

Aeq(2:nu+1,nx+nu+2:nx+nu+1+ncx+ncu) = [B'*H' Hu']; % Derivative of Lagrangian w.r.t. u = zero
%RHS already initialized to zero

%Objective function
cvec = [zeros(1,nx+nu) -1 zeros(1,2*(ncx+ncu))]';

intc=[repmat('C',nx+nu+1+ncx+ncu,1);repmat('B',ncx+ncu,1)]'; % C for continuous variables, B for binary

opts = cplexoptimset('mip.tolerances.mipgap',0.01,'mip.tolerances.absmipgap',0.01);%Default 1e-4,1e-6

[x,lam1,flag,ef] = cplexmilp(cvec,Aineq,bineq,Aeq,beq,[],[],[],blx,bux,intc,[],opts);




