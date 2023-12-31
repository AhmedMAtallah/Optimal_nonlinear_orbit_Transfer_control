
%     This file is part of CasADi.
%
%     CasADi -- A symbolic framework for dynamic optimization.
%     Copyright (C) 2010-2014 Joel Andersson, Joris Gillis, Moritz Diehl,
%                             K.U. Leuven. All rights reserved.
%     Copyright (C) 2011-2014 Greg Horn
%
%     CasADi is free software; you can redistribute it and/or
%     modify it under the terms of the GNU Lesser General Public
%     License as published by the Free Software Foundation; either
%     version 3 of the License, or (at your option) any later version.
%
%     CasADi is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
%     Lesser General Public License for more details.
%
%     You should have received a copy of the GNU Lesser General Public
%     License along with CasADi; if not, write to the Free Software
%     Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%
%
clear all
import casadi.*
% Declare variables (use simple, efficient DAG)
x0=SX.sym('x0'); x1=SX.sym('x1');
x = [x0;x1];

% Control
u = SX.sym('u');

% ODE right hand side
xdot = vertcat((1 - x1*x1)*x0 - x1 + u, x0);

% Lagrangian function
L = x0*x0 + x1*x1 + u*u;

% Costate
lam = SX.sym('lam',2);

% Hamiltonian function
H = dot(lam,xdot) + L;

% Costate equations
ldot = -gradient(H,x);

%% The control must minimize the Hamiltonian, which is:
disp(['Hamiltonian: ',str(H)]);


% H is of a convex quadratic form in u: H = u*u + p*u + q, let's get the coefficient p
p = gradient(H,u);     % this gives us 2*u + p
p = substitute(p,u,0); % replace u with zero: gives us p

% H's unconstrained minimizer is: u = -p/2
u_opt = -p/2;

% We must constrain u to the interval [-0.75, 1.0], convexity of H ensures that the optimum is obtain at the bound when u_opt is outside the interval
u_opt = fmin(u_opt,1.0);
u_opt = fmax(u_opt,-0.75);
disp(['optimal control: ',str(u_opt)])

% Augment f with lam_dot and subtitute in the value for the optimal control
f = [xdot;ldot];
f = substitute(f,u,u_opt);

% Function for obtaining the optimal control from the augmented state
u_fcn = Function('ufcn', {[x;lam]}, {u_opt});

% Formulate the DAE
dae = struct('x',[x;lam], 'ode',f);

% Augmented DAE state dimension
nX = 4;

% End time
tf = 10.0;

% Number of shooting nodes
num_nodes = 20;

% Create an integrator (CVodes)
iopts = struct;
iopts.abstol = 1e-8; % abs. tolerance
iopts.reltol = 1e-8; % rel. tolerance
iopts.t0 = 0.0;
iopts.tf = tf/num_nodes;
I = integrator('I', 'cvodes', dae, iopts);

% Variables for the states at each shooting node
X = MX.sym('X',nX,num_nodes+1);

% Formulate the root finding problem
G = {};
G{end+1} = X(1:2,1) - [0;1]; % states fixed, costates free at initial time
for k = 1:num_nodes
  XF = getfield(I('x0',X(:,k)),'xf');
  G{end+1} = XF-X(:,k+1);
end
G{end+1} = X(3:4,num_nodes+1) - [0;0]; % costates fixed, states free at final time

% Terminal constraints: lam = 0
rfp = struct('x',vec(X),'g',vertcat(G{:}));

% Select a solver for the root-finding problem
Solver = 'nlpsol';
%Solver = "newton"
%Solver = "kinsol"

% Solver options
opts = struct;
if Solver=='nlpsol',
    opts.nlpsol = 'ipopt';
    opts.nlpsol_options.ipopt.hessian_approximation = 'limited-memory';
elseif Solver=='newton',
    opts.linear_solver = 'csparse';
elseif Solver=='kinsol',
    opts.linear_solver_type = 'user_defined';
    opts.linear_solver = 'csparse';
    opts.max_iter = 1000;
end

% Allocate a solver
solver = rootfinder('solver', Solver, rfp, opts);

% Solve the problem
X_sol = getfield(solver('x0',0),'x');

% Time grid for visualization
tgrid = linspace(0,tf,100);

% Simulator to get optimal state and control trajectories
simulator = integrator('simulator', 'cvodes', dae, struct('grid',tgrid,'output_t0',true));

% Simulate to get the trajectories
sol = full(getfield(simulator('x0',X_sol(1:4)),'xf'));

% Calculate the optimal control
u_opt = full(u_fcn(sol));

% Plot the results
figure();
hold on
plot(tgrid, sol(1,:)', '--');
plot(tgrid, sol(2,:)', '-');
plot(tgrid, u_opt', '-.');
title('Van der Pol optimization - indirect multiple shooting');
xlabel('time');
legend('x trajectory','y trajectory','u trajectory');
grid on;