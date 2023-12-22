% Define the problem parameters
mu = 3.986e5;          % Earth's gravitational parameter (km^3/s^2)
R = 6378;              % Earth's radius (km)
m0 = 1000;             % Initial mass of spacecraft (kg)
mf = 800;              % Final mass of spacecraft (kg)
g0 = 9.81;             % Acceleration due to gravity at Earth's surface (m/s^2)
Isp = 300;             % Specific impulse of thruster (s)
Tmax = 2*Isp*g0*m0;    % Maximum thrust magnitude (N)

% Define the final state
rf = [2*R, 0, 0];      % Final position (km)
vf = [0, sqrt(mu/norm(rf)), 0]; % Final velocity (km/s)

% Define the initial guess for the control input
N = 100;               % Number of control intervals
t = linspace(0, 1000, N+1);     % Time vector (s)
u_guess = zeros(3,N);  % Initial guess for control input (N)

% Define the anonymous functions for the state derivative and the constraint
f = @(t,x,u) [x(4:6); -mu*x(1:3)/norm(x(1:3))^3 + u/m0];
g = @(t,x,u) 1/(g0*Isp)*norm(u);

% Define the boundary conditions
x0 = [R, 0, 0, 0, sqrt(mu/R), m0];  % Initial state
xf = [rf, vf, mf];                   % Final state

% Define the shooting function
shooting_fun = @(u_guess) orbit_transfer_shooting_fun(u_guess, t, f, g, x0, xf, Tmax);

% Use fsolve to find the optimal control input
options = optimoptions('fsolve', 'Display', 'iter');
u_opt = fsolve(shooting_fun, u_guess, options);

% Integrate the equations of motion and the costate equations using ode45
[t, x] = ode45(@(t,x) f(t,x,interp1(t,u_opt',t)'), [0, t(end)], x0);
[t, lambda] = ode45(@(t,l) costate_eqn(t,l,interp1(t,x(:,7)',t)',interp1(t,u_opt',t)',g0,Isp), [t(end), 0], [0,0,0,0,0,1]);

% Compute the optimal cost and state/costate vectors
J_opt = t(end);
x_opt = [x(:,1:6), lambda(:,1:6)];
u_opt = u_opt';

% Plot the results
figure;
subplot(4,2,1); plot(t, x(:,1)); title('Position x'); ylabel('km');
subplot(4,2,2); plot(t, x(:,2)); title('Position y'); ylabel('km');
subplot(4,2,3); plot(t, x(:,3)); title('Position z'); ylabel('km');
subplot(4,2,4); plot(t, x(:,4)); title('Velocity x'); ylabel('km/s');
subplot(4,2,5); plot(t, x(:,5)); title('Velocity y'); ylabel('km/s');
subplot(4,2,6); plot(t, x(:,6)); title('Velocity z')
